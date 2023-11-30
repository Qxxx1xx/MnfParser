#pragma once

#include <Eigen/Dense>
#include <bit>
#include <fstream>
#include <string>

namespace MnfParser {
template<typename T>
T
read_big_endian_value(std::istream& is)
{
    // 声明一个类型为 T 的变量 value，用于存储从文件中读取的数据
    T value;
    // 将 value 的地址转换为 char* 类型，read函数需要一个 char* 来读取数据
    char* buffer = reinterpret_cast<char*>(&value);
    // 从文件中读取 sizeof(T) 个字节，存储到 value 的地址中
    if (!is.read(buffer, sizeof(T))) {
        throw std::runtime_error("read_big_endian_value failed");
    }
    // 如果原生字节序是小端的，这行代码将反转从文件读取的字节
    if constexpr (std::endian::native == std::endian::little) {
        std::reverse(buffer, buffer + sizeof(T));
    }
    return value;
}

struct MnfBinParser
{
    std::ifstream mnf_file;

    // 常量数据成员
    // 整型数据类型的字节数
    const int kIntByteNum = 4;
    const int kDoubleByteNum = 8;
    // 版本相关
    const int kVersionByteOffset = 0x3;
    const int kSupportMnfVersion = 6;
    const int kVersionMultiplier = 10;
    // 节点坐标相关
    const int kNodeNumByteOffset = 0xa10 + 4;
    const int kCoordDimByteOffset = kNodeNumByteOffset + kIntByteNum;
    const int kCoordDim3NodeDof = 6;

    enum class EigenValuesType
    {
        kEigenValue,
        kNaturalFreqHz,
        kNaturalFreqRad,
        kGeneralizedMass,
        kGeneralizedStiffness,
    };

    explicit MnfBinParser(const std::string& file_path);

    void checkMnfVersion()
    {
        char char_mnf_version;
        int expected_mnf_version = kSupportMnfVersion * kVersionMultiplier;

        mnf_file.seekg(kVersionByteOffset, std::ios::beg);
        mnf_file.read(&char_mnf_version, 1);

        if (!mnf_file) {
            throw std::runtime_error("Failed to read MNF version from file");
        }

        if (char_mnf_version != expected_mnf_version) {
            std::string errMsg = "Invalid MNF Version (" + std::to_string(char_mnf_version / 10.0) +
                                 ") - This solver only supports MNF Version " + std::to_string(kSupportMnfVersion);
            throw std::runtime_error(errMsg);
        }
    }

    int getNodeNum()
    {
        mnf_file.seekg(kNodeNumByteOffset, std::ios::beg);
        return read_big_endian_value<int>(mnf_file);
    }

    int getCoordDim()
    {
        mnf_file.seekg(kCoordDimByteOffset, std::ios::beg);
        return read_big_endian_value<int>(mnf_file);
    }

    int byteOffsetToNodeCoord() { return kCoordDimByteOffset + kIntByteNum; }

    int byteNumToNodeCoord()
    {
        int node_num = getNodeNum();
        int coord_dim = getCoordDim();
        return node_num * (kIntByteNum + coord_dim * sizeof(double));
    }

    Eigen::MatrixXd getNodeCoord()
    {
        int node_num = getNodeNum();
        int coord_dim = getCoordDim();
        Eigen::MatrixXd node_coord(node_num, coord_dim);

        // 节点坐标相关
        const int kNodeCoordByteOffset = byteOffsetToNodeCoord();

        mnf_file.seekg(kNodeCoordByteOffset, std::ios::beg);
        for (int i = 0; i < node_num; ++i) {
            // 跳过节点序号的整型数据类型的字节数
            mnf_file.seekg(kIntByteNum, std::ios::cur);
            for (int j = 0; j < coord_dim; ++j) {
                node_coord(i, j) = read_big_endian_value<double>(mnf_file);
            }
        }

        if (!mnf_file) {
            throw std::runtime_error("Failed to read node coordinates from file");
        }

        return node_coord;
    }

    int byteOffsetToMass() { return byteOffsetToNodeCoord() + byteNumToNodeCoord(); }

    int byteOffsetToCenterOfMass() { return byteOffsetToMass() + kDoubleByteNum; }

    constexpr int byteNumToCenterOfMass() { return kDoubleByteNum * 3; }

    int byteOffsetToInertiaTensor() { return byteOffsetToCenterOfMass() + byteNumToCenterOfMass(); }

    constexpr int byteNumToInertiaTensor() { return kDoubleByteNum * 9; }

    double getMass()
    {
        mnf_file.seekg(byteOffsetToMass(), std::ios::beg);
        return read_big_endian_value<double>(mnf_file);
    }

    Eigen::Vector3d getCenterOfMass()
    {
        mnf_file.seekg(byteOffsetToCenterOfMass(), std::ios::beg);
        Eigen::Vector3d center_of_mass;
        for (int i = 0; i < 3; ++i) {
            center_of_mass(i) = read_big_endian_value<double>(mnf_file);
        }
        return center_of_mass;
    }

    Eigen::Matrix3d getInertiaTensor()
    {
        mnf_file.seekg(byteOffsetToInertiaTensor(), std::ios::beg);
        Eigen::Matrix3d inertia_tensor;
        for (int i = 0; i < 9; ++i) {
            inertia_tensor(i) = read_big_endian_value<double>(mnf_file);
        }
        return inertia_tensor;
    }

    int byteOffsetToModalOrder() { return byteOffsetToInertiaTensor() + byteNumToInertiaTensor(); }

    int getModalOrder()
    {
        mnf_file.seekg(byteOffsetToModalOrder(), std::ios::beg);
        return read_big_endian_value<int>(mnf_file);
    }

    int byteOffsetToEigenvalues() { return byteOffsetToModalOrder() + kIntByteNum; }

    Eigen::VectorXd _getEigenvalues(EigenValuesType eigen_values_type)
    {
        int modal_order = getModalOrder();
        Eigen::VectorXd eigenvalues(modal_order);

        // 得到eigen_values_type在EigenValuesType位置
        int eigen_values_type_position = static_cast<int>(eigen_values_type);

        int pre_byte_offest = kIntByteNum + eigen_values_type_position * kDoubleByteNum;
        int rear_byte_offest = (5 - eigen_values_type_position - 1) * kDoubleByteNum;

        mnf_file.seekg(byteOffsetToEigenvalues(), std::ios::beg);
        for (int i = 0; i < modal_order; ++i) {
            // 跳过前面的字节
            mnf_file.seekg(pre_byte_offest, std::ios::cur);
            eigenvalues(i) = read_big_endian_value<double>(mnf_file);
            // 跳过后面的字节
            mnf_file.seekg(rear_byte_offest, std::ios::cur);
        }
        return eigenvalues;
    }

    Eigen::VectorXd getEigenvalues() { return _getEigenvalues(EigenValuesType::kEigenValue); }

    Eigen::VectorXd getNaturalFreqHz() { return _getEigenvalues(EigenValuesType::kNaturalFreqHz); }

    Eigen::VectorXd getNaturalFreqRad() { return _getEigenvalues(EigenValuesType::kNaturalFreqRad); }

    Eigen::VectorXd getGeneralizedMass() { return _getEigenvalues(EigenValuesType::kGeneralizedMass); }

    Eigen::VectorXd getGeneralizedStiffness() { return _getEigenvalues(EigenValuesType::kGeneralizedStiffness); }

    // 计算节点自由度
    int getNodeDof()
    {
        int coord_dim = getCoordDim();
        if (coord_dim == 3) {
            return kCoordDim3NodeDof;
        }
        // todo: 未实现其他维度的节点自由度
        return 0;
    }

    int byteNumToEigenvalues()
    {
        int modal_order = getModalOrder();
        return (kIntByteNum + 5 * kDoubleByteNum) * modal_order;
    }

    int byteOffsetToModeShapes()
    {
        int modal_order = getModalOrder();
        int node_num = getNodeNum();
        return byteOffsetToEigenvalues() + byteNumToEigenvalues();
    }

    Eigen::MatrixXd getModeShapes()
    {
        int node_num = getNodeNum();
        int coord_dim = getCoordDim();
        int modal_order = getModalOrder();
        // 节点总自由度
        int nodes_total_dof = node_num * getNodeDof();
        Eigen::MatrixXd mode_shapes(nodes_total_dof, modal_order);

        mnf_file.seekg(byteOffsetToModeShapes(), std::ios::beg);
        if (read_big_endian_value<int>(mnf_file) != modal_order) {
            throw std::runtime_error("Failed to read mode shapes from file");
        }
        if (read_big_endian_value<int>(mnf_file) != node_num) {
            throw std::runtime_error("Failed to read mode shapes from file");
        }
        int node_num_temp = 0;
        for (int i = 0; i < node_num; i++) {
            node_num_temp = read_big_endian_value<int>(mnf_file);
        }
        if (node_num_temp != node_num) {
            throw std::runtime_error("Failed to read mode shapes from file");
        }

        for (int i = 0; i < modal_order; i++) {
            mnf_file.seekg(kIntByteNum, std::ios::cur);
            for (int j = 0; j < nodes_total_dof; j++) {
                mode_shapes(j, i) = read_big_endian_value<double>(mnf_file);
            }
        }

        if (!mnf_file) {
            throw std::runtime_error("Failed to read mode shapes from file");
        }

        return mode_shapes;
    }

    int byteNumToModeShapes()
    {
        int node_num = getNodeNum();
        int modal_order = getModalOrder();
        int nodes_total_dof = node_num * getNodeDof();
        return kIntByteNum + kIntByteNum + kIntByteNum * node_num + modal_order * kIntByteNum +
               modal_order * nodes_total_dof * kDoubleByteNum;
    }

    int byteOffsetToNodalMasses() { return byteOffsetToModeShapes() + byteNumToModeShapes(); }

    Eigen::VectorXd getNodalMasses()
    {
        int node_num = getNodeNum();
        Eigen::VectorXd nodal_masses(node_num);

        mnf_file.seekg(byteOffsetToNodalMasses(), std::ios::beg);
        if (read_big_endian_value<int>(mnf_file) != node_num) {
            throw std::runtime_error("Failed to read nodal masses from file");
        }
        if (read_big_endian_value<int>(mnf_file) != 1) {
            throw std::runtime_error("Failed to read nodal masses from file");
        }

        for (int i = 0; i < node_num; ++i) {
            mnf_file.seekg(kIntByteNum, std::ios::cur);
            nodal_masses(i) = read_big_endian_value<double>(mnf_file);
        }

        if (!mnf_file) {
            throw std::runtime_error("Failed to read nodal masses from file");
        }

        return nodal_masses;
    }

    int byteNumToNodalMasses() { return kIntByteNum + kIntByteNum + getNodeNum() * (kIntByteNum + kDoubleByteNum); }

    int byteOffsetToNodalInertias() { return byteOffsetToNodalMasses() + byteNumToNodalMasses(); }

    std::vector<Eigen::Vector3d> getNodalInertias()
    {
        int node_num = getNodeNum();
        std::vector<Eigen::Vector3d> nodal_inertias(node_num);

        mnf_file.seekg(byteOffsetToNodalInertias(), std::ios::beg);
        int entry_num = read_big_endian_value<int>(mnf_file);

        for (int i = 0; i < node_num; ++i) {
            for (int j = 0; j < 3; ++j) {
                mnf_file.seekg(kIntByteNum + kIntByteNum + kIntByteNum, std::ios::cur);
                nodal_inertias[i](j) = read_big_endian_value<double>(mnf_file);
            }
        }

        if (!mnf_file) {
            throw std::runtime_error("Failed to read nodal inertias from file");
        }

        return nodal_inertias;
    }
};
} // namespace MnfParser