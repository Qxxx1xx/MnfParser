#include <MnfParser/MnfBinParser.h>
#include <gtest/gtest.h>
#include <string>

struct MnfBinParserTest : public ::testing::Test
{
    MnfParser::MnfBinParser mnf_bin_parser;
    const double kTolerance = 1e-9;

    MnfBinParserTest()
        : mnf_bin_parser(std::string(MnfParserTESTS_DATA_DIR) + "siso_beam.mnf")
    {
        // 这里可以添加更多初始化代码
    }

    virtual ~MnfBinParserTest()
    {
        // 这里可以添加清理代码
    }
};

TEST_F(MnfBinParserTest, Test_MnfParserTESTS_DATA_DIR)
{
    EXPECT_EQ(std::string(MnfParserTESTS_DATA_DIR), "C:/Users/Qxxxx/Desktop/Work/code/MnfParser/tests/data/");
}

// MnfBinParser构造函数测试
TEST_F(MnfBinParserTest, MnfBinParserConstructor)
{
    EXPECT_EQ(mnf_bin_parser.mnf_file.is_open(), true);
}

// MnfBinParser::getNodeNum()测试
TEST_F(MnfBinParserTest, getNodeNum)
{
    EXPECT_EQ(mnf_bin_parser.getNodeNum(), 101);
}

// MnfBinParser::getCoordDim()测试
TEST_F(MnfBinParserTest, getCoordDim)
{
    EXPECT_EQ(mnf_bin_parser.getCoordDim(), 3);
}

// MnfBinParser::getNodeCoord()测试
TEST_F(MnfBinParserTest, getNodeCoord)
{
    Eigen::MatrixXd node_coord = mnf_bin_parser.getNodeCoord();
    EXPECT_EQ(node_coord.rows(), 101);
    EXPECT_EQ(node_coord.cols(), 3);
    EXPECT_EQ(double(node_coord(1, 0)), 1.0);
}

// MnfBinParser::getMass()测试
TEST_F(MnfBinParserTest, getMass)
{
    double mass = mnf_bin_parser.getMass();
    EXPECT_NEAR(mass, 0.27667, kTolerance);
}

// MnfBinParser::getCenterOfMass()测试
TEST_F(MnfBinParserTest, getCenterOfMass)
{
    Eigen::Vector3d center_of_mass = mnf_bin_parser.getCenterOfMass();
    EXPECT_NEAR(center_of_mass(0), 0.5, kTolerance);
    EXPECT_NEAR(center_of_mass(1), 0.0, kTolerance);
    EXPECT_NEAR(center_of_mass(2), 0.0, kTolerance);
}

// MnfBinParser::getInertiaTensor()测试
TEST_F(MnfBinParserTest, getInertiaTensor)
{
    Eigen::Matrix3d inertia_tensor = mnf_bin_parser.getInertiaTensor();
    EXPECT_NEAR(inertia_tensor(0, 0), 0.0000046112, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 1), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(0, 2), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 0), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 1), 0.0922325557, kTolerance);
    EXPECT_NEAR(inertia_tensor(1, 2), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 0), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 1), 0.0, kTolerance);
    EXPECT_NEAR(inertia_tensor(2, 2), 0.0922325557, kTolerance);
}

// MnfBinParser::_getEigenVectors()测试
TEST_F(MnfBinParserTest, _getEigenVectors)
{
    int modal_order = mnf_bin_parser.getModalOrder();
    EXPECT_EQ(modal_order, 18);
    Eigen::VectorXd eigenvectors = mnf_bin_parser.getEigenvalues();
    EXPECT_NEAR(eigenvectors(0), 3.7869932612e-07, kTolerance);
    Eigen::VectorXd freq_hz = mnf_bin_parser.getNaturalFreqHz();
    EXPECT_NEAR(freq_hz(0), 9.7941645850e-05, kTolerance);
    Eigen::VectorXd freq_rad = mnf_bin_parser.getNaturalFreqRad();
    EXPECT_NEAR(freq_rad(0), 6.1538551016e-04, kTolerance);
    Eigen::VectorXd generalized_mass = mnf_bin_parser.getGeneralizedMass();
    EXPECT_NEAR(generalized_mass(0), 1.0, kTolerance);
    Eigen::VectorXd generalized_stiffness = mnf_bin_parser.getGeneralizedStiffness();
    EXPECT_NEAR(generalized_stiffness(0), 3.7869932612e-07, kTolerance);
}

// MnfBinParser::getModeShapes()测试
TEST_F(MnfBinParserTest, getModeShapes)
{
    Eigen::MatrixXd mode_shapes = mnf_bin_parser.getModeShapes();
    EXPECT_NEAR(mode_shapes(0, 17), -3.292578719e+00, kTolerance);
}

// MnfBinParser::getNodalMasses()测试
TEST_F(MnfBinParserTest, getNodalMasses)
{
    Eigen::VectorXd nodal_masses = mnf_bin_parser.getNodalMasses();
    EXPECT_NEAR(nodal_masses(0), 1.383350e-03, kTolerance);
    EXPECT_NEAR(nodal_masses(100), 2.766700e-03, kTolerance);
}

// MnfBinParser::getNodalInertias()测试
TEST_F(MnfBinParserTest, getNodalInertias)
{
    std::vector<Eigen::Vector3d> nodal_inertias = mnf_bin_parser.getNodalInertias();
    EXPECT_NEAR(nodal_inertias[0](0), 0.0000000231, kTolerance);
    EXPECT_NEAR(nodal_inertias[2](0), 0.0000000461, kTolerance);
}

// MnfBinParser::getElementFaces()测试
TEST_F(MnfBinParserTest, getElementFaces)
{
    int faces_num = mnf_bin_parser.getFacesNum();
    EXPECT_EQ(faces_num, 100);
    int faces_data_int_num = mnf_bin_parser.getFacesDataIntNum();
    EXPECT_EQ(faces_data_int_num, 300);
    auto element_faces = mnf_bin_parser.getElementFaces();
    EXPECT_EQ(element_faces[0][0], 1);
    EXPECT_EQ(element_faces[0][1], 3);
    EXPECT_EQ(element_faces.back()[0], 101);
    EXPECT_EQ(element_faces.back()[1], 2);
}

// MnfBinParser::getModeShapeTransformation()测试
TEST_F(MnfBinParserTest, getModeShapeTransformation)
{
    Eigen::MatrixXd mode_shape_transformation = mnf_bin_parser.getModeShapeTransformation();
    int modal_order = mnf_bin_parser.getModalOrder();
    // 判断是否是modal_order * modal_order的单位矩阵
    bool is_identity_matrix = mode_shape_transformation.isIdentity();
    EXPECT_EQ(is_identity_matrix, true);
    EXPECT_EQ(mode_shape_transformation.rows(), modal_order);
}

// MnfBinParser::getInterfaceNodes()测试
TEST_F(MnfBinParserTest, getInterfaceNodes)
{
    int interface_nodes_num = mnf_bin_parser.getInterfaceNodesNum();
    EXPECT_EQ(interface_nodes_num, 2);
    auto interface_nodes = mnf_bin_parser.getInterfaceNodes();
    EXPECT_EQ(interface_nodes[0], 1);
    EXPECT_EQ(interface_nodes[1], 2);
}

// mnf_bin_parser.mnf_file.seekg(mnf_bin_parser.byteOffsetToElementFaces(), std::ios::beg);
// std::cout << MnfParser::read_big_endian_value<int>(mnf_bin_parser.mnf_file) << std::endl;
