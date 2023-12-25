#include <MnfParser/MnfBinParser.h>
#include <gtest/gtest.h>
#include <string>

struct MnfBinParserTestSisoBlock : public ::testing::Test
{
    MnfParser::MnfBinParser mnf_bin_parser;
    const double kTolerance = 1e-9;

    MnfBinParserTestSisoBlock()
        : mnf_bin_parser(std::string(MnfParserTESTS_DATA_DIR) + "siso_block.mnf")
    {
        // 这里可以添加更多初始化代码
    }

    virtual ~MnfBinParserTestSisoBlock()
    {
        // 这里可以添加清理代码
    }
};

TEST_F(MnfBinParserTestSisoBlock, Test_MnfParserTESTS_DATA_DIR)
{
    EXPECT_EQ(std::string(MnfParserTESTS_DATA_DIR), "E:/Work/code/MnfParser/tests/data/");
}

// MnfBinParser构造函数测试
TEST_F(MnfBinParserTestSisoBlock, MnfBinParserConstructor)
{
    EXPECT_EQ(mnf_bin_parser.mnf_file.is_open(), true);
}

// MnfBinParser::hasNodalInertias()测试
TEST_F(MnfBinParserTestSisoBlock, hasNodalInertias)
{
    EXPECT_EQ(mnf_bin_parser.hasNodalInertias(), true);
}

// MnfBinParser::getNodeNum()测试
TEST_F(MnfBinParserTestSisoBlock, getNodeNum)
{
    EXPECT_EQ(mnf_bin_parser.getNodeNum(), 191);
}

// MnfBinParser::getCoordDim()测试
TEST_F(MnfBinParserTestSisoBlock, getCoordDim)
{
    EXPECT_EQ(mnf_bin_parser.getCoordDim(), 3);
}

// MnfBinParser::getNodeCoord()测试
TEST_F(MnfBinParserTestSisoBlock, getNodeCoord)
{
    Eigen::MatrixXd node_coord = mnf_bin_parser.getNodeCoord();
    EXPECT_EQ(node_coord.rows(), 191);
    EXPECT_EQ(node_coord.cols(), 3);
    EXPECT_NEAR(double(node_coord(190, 0)), 5.0000000000e-02, kTolerance);
}

// MnfBinParser::getMass()测试
TEST_F(MnfBinParserTestSisoBlock, getMass)
{
    double mass = mnf_bin_parser.getMass();
    EXPECT_NEAR(mass, 27.0, kTolerance);
}

// MnfBinParser::getCenterOfMass()测试
TEST_F(MnfBinParserTestSisoBlock, getCenterOfMass)
{
    Eigen::Vector3d center_of_mass = mnf_bin_parser.getCenterOfMass();
    EXPECT_NEAR(center_of_mass(0), 0.05, kTolerance);
    EXPECT_NEAR(center_of_mass(1), 0.05, kTolerance);
    EXPECT_NEAR(center_of_mass(2), 0.5, kTolerance);
}

// MnfBinParser::getInertiaTensor()测试
TEST_F(MnfBinParserTestSisoBlock, getInertiaTensor)
{
    Eigen::Matrix3d inertia_tensor = mnf_bin_parser.getInertiaTensor();
    EXPECT_NEAR(inertia_tensor(0, 0), 9.1125, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 1), 0.0675, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 2), 0.6750, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 0), 0.0675, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 1), 9.1125, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 2), 0.6750, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 0), 0.6750, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 1), 0.6750, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 2), 0.2025, kTolerance);
}

// MnfBinParser::_getEigenVectors()测试
TEST_F(MnfBinParserTestSisoBlock, _getEigenVectors)
{
    int modal_order = mnf_bin_parser.getModalOrder();
    EXPECT_EQ(modal_order, 18);
    Eigen::VectorXd eigenvectors = mnf_bin_parser.getEigenvalues();
    EXPECT_NEAR(eigenvectors(6), 7.1478616400e+06, 1.0);
    Eigen::VectorXd freq_hz = mnf_bin_parser.getNaturalFreqHz();
    EXPECT_NEAR(freq_hz(6), 4.2550846110e+02, 1.0);
    Eigen::VectorXd freq_rad = mnf_bin_parser.getNaturalFreqRad();
    EXPECT_NEAR(freq_rad(6), 2.6735485109e+03, 1.0);
    Eigen::VectorXd generalized_mass = mnf_bin_parser.getGeneralizedMass();
    EXPECT_NEAR(generalized_mass(6), 1.0, kTolerance);
    Eigen::VectorXd generalized_stiffness = mnf_bin_parser.getGeneralizedStiffness();
    EXPECT_NEAR(generalized_stiffness(6), 7.1478616400e+06, 1.0);
}

// MnfBinParser::getModeShapes()测试
TEST_F(MnfBinParserTestSisoBlock, getModeShapes)
{
    Eigen::MatrixXd mode_shapes = mnf_bin_parser.getModeShapes();
    EXPECT_NEAR(mode_shapes(0, 17), 2.600650995e-01, kTolerance);
}

// MnfBinParser::getNodalMasses()测试
TEST_F(MnfBinParserTestSisoBlock, getNodalMasses)
{
    Eigen::VectorXd nodal_masses = mnf_bin_parser.getNodalMasses();
    EXPECT_NEAR(nodal_masses(0), 4.218750e-02, kTolerance);
    EXPECT_NEAR(nodal_masses(100), 8.437500e-02, kTolerance);
}

// MnfBinParser::getNodesAndNodalInertias()测试
TEST_F(MnfBinParserTestSisoBlock, getNodesAndNodalInertias)
{
    auto [nodes, nodal_inertias] = mnf_bin_parser.getNodesAndNodalInertias();
    EXPECT_EQ(nodes.size(), 2);
    EXPECT_EQ(nodal_inertias.size(), 2);
    EXPECT_NEAR(nodal_inertias[0](0), 0.0, kTolerance);
    EXPECT_NEAR(nodal_inertias[1](0), 0.0, kTolerance);
}

// MnfBinParser::getElementFaces()测试
TEST_F(MnfBinParserTestSisoBlock, getElementFaces)
{
    int faces_num = mnf_bin_parser.getFacesNum();
    EXPECT_EQ(faces_num, 480);
    int faces_data_int_num = mnf_bin_parser.getFacesDataIntNum();
    EXPECT_EQ(faces_data_int_num, 480 * 5);
    auto element_faces = mnf_bin_parser.getElementFaces();
    EXPECT_EQ(element_faces[0][0], 2);
    EXPECT_EQ(element_faces[0][1], 3);
    EXPECT_EQ(element_faces.back()[0], 18);
    EXPECT_EQ(element_faces.back()[1], 14);
}

// MnfBinParser::getModeShapeTransformation()测试
TEST_F(MnfBinParserTestSisoBlock, getModeShapeTransformation)
{
    Eigen::MatrixXd mode_shape_transformation = mnf_bin_parser.getModeShapeTransformation();
    int modal_order = mnf_bin_parser.getModalOrder();
    // 判断是否是modal_order * modal_order的单位矩阵
    bool is_identity_matrix = mode_shape_transformation.isIdentity();
    EXPECT_EQ(is_identity_matrix, true);
    EXPECT_EQ(mode_shape_transformation.rows(), modal_order);
}

// MnfBinParser::getInterfaceNodes()测试
TEST_F(MnfBinParserTestSisoBlock, getInterfaceNodes)
{
    int interface_nodes_num = mnf_bin_parser.getInterfaceNodesNum();
    EXPECT_EQ(interface_nodes_num, 2);
    auto interface_nodes = mnf_bin_parser.getInterfaceNodes();
    EXPECT_EQ(interface_nodes[0], 190);
    EXPECT_EQ(interface_nodes[1], 191);
}

// mnf_bin_parser.mnf_file.seekg(mnf_bin_parser.byteOffsetToElementFaces(), std::ios::beg);
// std::cout << MnfParser::read_big_endian_value<int>(mnf_bin_parser.mnf_file) << std::endl;
