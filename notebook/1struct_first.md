设计一个类来存储和处理这些数据，需要考虑数据的封装、类的职责划分以及易于维护和扩展的设计原则。下面是一个基本的设计草案，用于处理你提供的功能：

### 类设计概览

1. **类名称**：以功能命名，如 `MnfParser` 或 `StructuralDataLoader`。
2. **私有成员变量**：存储文件路径、节点数据、质量数据、模态数据等。
3. **公共方法**：包括构造函数、文件读取函数、数据访问函数等。
4. **异常处理**：在类内部处理文件读取和数据处理中的异常。

### 类设计细节

#### 1. 类定义

```cpp
class StructuralDataLoader {
public:
    StructuralDataLoader(const std::string& filePath);
    void readMnfFile();

    // 数据访问方法
    std::vector<Node> getNodes() const;
    // ... 其他数据访问方法

private:
    std::string filePath;
    std::vector<Node> nodes;
    std::vector<ModalData> modalData;
    // ... 其他必要的成员变量

    // 读取和解析文件的辅助方法
    void parseVersionInfo(std::ifstream& file);
    void parseNodeCoords(std::ifstream& file);
    void parseGlobalProperties(std::ifstream& file);
    // ... 其他解析方法
};
```

#### 2. 构造函数和文件读取

```cpp
StructuralDataLoader::StructuralDataLoader(const std::string& filePath) : filePath(filePath) {
}

void StructuralDataLoader::readMnfFile() {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath);
    }

    parseVersionInfo(file);
    parseNodeCoords(file);
    parseGlobalProperties(file);
    // ... 调用其他解析方法
}
```

#### 3. 解析方法

```cpp
void StructuralDataLoader::parseVersionInfo(std::ifstream& file) {
    // 实现版本信息的解析
}

void StructuralDataLoader::parseNodeCoords(std::ifstream& file) {
    // 实现节点坐标的解析
}

void StructuralDataLoader::parseGlobalProperties(std::ifstream& file) {
    // 实现全局属性的解析
}
```

#### 4. 数据访问方法

```cpp
std::vector<Node> StructuralDataLoader::getNodes() const {
    return nodes;
}

// ... 其他数据访问方法
```

### 注意事项

- **封装**：类应该封装所有与读取和解析文件相关的逻辑，只暴露必要的接口给用户。
- **异常处理**：确保在读取文件或解析数据时遇到错误能够抛出并处理异常。
- **可扩展性**：设计时考虑到未来可能对类进行扩展，比如添加新的数据解析方法或改进现有方法。
- **代码组织**：相关的功能应该组织在一起，比如所有解析相关的方法放在一起。
- **文档和注释**：提供清晰的文档和注释，说明每个方法的用途和工作方式。

这个设计是一个起点，根据实际需求和应用场景，可以进一步优化和调整。