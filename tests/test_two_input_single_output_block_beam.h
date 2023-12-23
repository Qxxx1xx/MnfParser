#include <MnfParser/MnfBinParser.h>
#include <gtest/gtest.h>
#include <string>

struct MnfBinParserTestMisoBlock : public ::testing::Test
{
    MnfParser::MnfBinParser mnf_bin_parser;
    const double kTolerance = 1e-9;

    MnfBinParserTestMisoBlock()
        : mnf_bin_parser(std::string(MnfParserTESTS_DATA_DIR) + "two_input_single_output_block_beam.mnf")
    {
        // 这里可以添加更多初始化代码
    }

    virtual ~MnfBinParserTestMisoBlock()
    {
        // 这里可以添加清理代码
    }
};

TEST_F(MnfBinParserTestMisoBlock, Test_MnfParserTESTS_DATA_DIR)
{
    EXPECT_EQ(std::string(MnfParserTESTS_DATA_DIR), "E:/Work/code/MnfParser/tests/data/");
}

// MnfBinParser构造函数测试
TEST_F(MnfBinParserTestMisoBlock, MnfBinParserConstructor)
{
    EXPECT_EQ(mnf_bin_parser.mnf_file.is_open(), true);
}

// MnfBinParser::hasNodalInertias()测试
TEST_F(MnfBinParserTestMisoBlock, hasNodalInertias)
{
    EXPECT_EQ(mnf_bin_parser.hasNodalInertias(), false);
}

// MnfBinParser::getNodeNum()测试
TEST_F(MnfBinParserTestMisoBlock, getNodeNum)
{
    EXPECT_EQ(mnf_bin_parser.getNodeNum(), 198);
}

// MnfBinParser::getCoordDim()测试
TEST_F(MnfBinParserTestMisoBlock, getCoordDim)
{
    EXPECT_EQ(mnf_bin_parser.getCoordDim(), 3);
}

// MnfBinParser::getNodeCoord()测试
TEST_F(MnfBinParserTestMisoBlock, getNodeCoord)
{
    Eigen::MatrixXd node_coord = mnf_bin_parser.getNodeCoord();
    EXPECT_EQ(node_coord.rows(), 198);
    EXPECT_EQ(node_coord.cols(), 3);
    EXPECT_NEAR(double(node_coord(1, 0)), 5.0000000000e-02, kTolerance);
}

// MnfBinParser::getMass()测试
TEST_F(MnfBinParserTestMisoBlock, getMass)
{
    double mass = mnf_bin_parser.getMass();
    EXPECT_NEAR(mass, 27.0000000000, kTolerance);
}

// MnfBinParser::getCenterOfMass()测试
TEST_F(MnfBinParserTestMisoBlock, getCenterOfMass)
{
    Eigen::Vector3d center_of_mass = mnf_bin_parser.getCenterOfMass();
    EXPECT_NEAR(center_of_mass(0), 0.0500000000, kTolerance);
    EXPECT_NEAR(center_of_mass(1), 0.0500000000, kTolerance);
    EXPECT_NEAR(center_of_mass(2), 0.5, kTolerance);
}

// MnfBinParser::getInertiaTensor()测试
TEST_F(MnfBinParserTestMisoBlock, getInertiaTensor)
{
    Eigen::Matrix3d inertia_tensor = mnf_bin_parser.getInertiaTensor();
    EXPECT_NEAR(inertia_tensor(0, 0), 9.1125000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 1), 0.0675000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 2), 0.6750000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 0), 0.0675000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 1), 9.1125000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 2), 0.6750000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 0), 0.6750000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 1), 0.6750000000, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 2), 0.2025000000, kTolerance);
}

// MnfBinParser::_getEigenVectors()测试
TEST_F(MnfBinParserTestMisoBlock, _getEigenVectors)
{
    int modal_order = mnf_bin_parser.getModalOrder();
    EXPECT_EQ(modal_order, 24);
    Eigen::VectorXd eigenvectors = mnf_bin_parser.getEigenvalues();
    EXPECT_NEAR(eigenvectors(9), 4.8093054462e+07, 1);
    Eigen::VectorXd freq_hz = mnf_bin_parser.getNaturalFreqHz();
    EXPECT_NEAR(freq_hz(9), 1.1037260986e+03, 1);
    Eigen::VectorXd freq_rad = mnf_bin_parser.getNaturalFreqRad();
    EXPECT_NEAR(freq_rad(9), 6.9349156060e+03, 1);
    Eigen::VectorXd generalized_mass = mnf_bin_parser.getGeneralizedMass();
    EXPECT_NEAR(generalized_mass(9), 1.0, kTolerance);
    Eigen::VectorXd generalized_stiffness = mnf_bin_parser.getGeneralizedStiffness();
    EXPECT_NEAR(generalized_stiffness(9), 4.8093054462e+07, 1);
}

// MnfBinParser::getModeShapes()测试
TEST_F(MnfBinParserTestMisoBlock, getModeShapes)
{
    Eigen::MatrixXd mode_shapes = mnf_bin_parser.getModeShapes();
    EXPECT_NEAR(mode_shapes(19 * 6, 17), -1.945906190e-02, kTolerance);
}

// MnfBinParser::getNodalMasses()测试
TEST_F(MnfBinParserTestMisoBlock, getNodalMasses)
{
    Eigen::VectorXd nodal_masses = mnf_bin_parser.getNodalMasses();
    EXPECT_NEAR(nodal_masses(0), 3.375000e-01, kTolerance);
    EXPECT_NEAR(nodal_masses(100), 4.218750e-02, kTolerance);
    EXPECT_NEAR(nodal_masses(197), 0.0, kTolerance);
}

// MnfBinParser::getNodesAndNodalInertias()测试
TEST_F(MnfBinParserTestMisoBlock, getNodesAndNodalInertias)
{
    auto [nodes, nodal_inertias] = mnf_bin_parser.getNodesAndNodalInertias();
    EXPECT_EQ(nodes.size(), 0);
    EXPECT_EQ(nodal_inertias.size(), 0);
}

// MnfBinParser::getElementFaces()测试
TEST_F(MnfBinParserTestMisoBlock, getElementFaces)
{
    int faces_num = mnf_bin_parser.getFacesNum();
    EXPECT_EQ(faces_num, 488);
    int faces_data_int_num = mnf_bin_parser.getFacesDataIntNum();
    EXPECT_EQ(faces_data_int_num, 488 * 5);
    auto element_faces = mnf_bin_parser.getElementFaces();
    EXPECT_EQ(element_faces[0][0], 1);
    EXPECT_EQ(element_faces[0][1], 188);
    EXPECT_EQ(element_faces.back()[0], 58);
    EXPECT_EQ(element_faces.back()[1], 81);
}

// MnfBinParser::getModeShapeTransformation()测试
TEST_F(MnfBinParserTestMisoBlock, getModeShapeTransformation)
{
    Eigen::MatrixXd mode_shape_transformation = mnf_bin_parser.getModeShapeTransformation();
    int modal_order = mnf_bin_parser.getModalOrder();
    // 判断是否是modal_order * modal_order的单位矩阵
    bool is_identity_matrix = mode_shape_transformation.isIdentity();
    EXPECT_EQ(is_identity_matrix, true);
    EXPECT_EQ(mode_shape_transformation.rows(), modal_order);
}

// MnfBinParser::getInterfaceNodes()测试
TEST_F(MnfBinParserTestMisoBlock, getInterfaceNodes)
{
    int interface_nodes_num = mnf_bin_parser.getInterfaceNodesNum();
    EXPECT_EQ(interface_nodes_num, 3);
    auto interface_nodes = mnf_bin_parser.getInterfaceNodes();
    EXPECT_EQ(interface_nodes[0], 190);
    EXPECT_EQ(interface_nodes[1], 191);
    EXPECT_EQ(interface_nodes[2], 192);
}

// mnf_bin_parser.mnf_file.seekg(mnf_bin_parser.byteOffsetToElementFaces(), std::ios::beg);
// std::cout << MnfParser::read_big_endian_value<int>(mnf_bin_parser.mnf_file) << std::endl;
