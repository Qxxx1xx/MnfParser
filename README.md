# MnfParser 项目

## 简介

MnfParser 是一个使用 C++20 构建的项目，目的是 [这里填写项目的简单描述和目的]。本项目使用 CMake 作为构建系统，并通过 vcpkg 管理依赖库。

## 开始

### 先决条件

为了构建和运行这个项目，你需要安装以下软件：

- CMake (3.15 或更高版本)
- 一个支持 C++20 的 C++ 编译器，如 GCC, Clang 或 MSVC
- Git，用于克隆仓库和下载 vcpkg

### 获取代码

首先，克隆项目到本地：

```bash
git clone https://github.com/Qxxx1xx/MnfParser.git
cd [项目目录]
```

### 使用 vcpkg 安装依赖

项目使用 vcpkg 来管理依赖。首次构建前，你需要安装 vcpkg 并安装必要的库：

```bash
./install_vcpkg.bat
```

这个脚本会自动从 GitHub 克隆 vcpkg 仓库并安装项目所需的所有依赖。

### 构建项目

使用 CMake 构建项目：

```bash
mkdir build
cmake -B build
cmake --build build
```

### 运行测试

项目支持自动测试，可以通过以下命令运行测试：

```bash
cd build
ctest -C Release
```

## 贡献

欢迎对 MnfParser 项目做出贡献。你可以通过 Pull Requests 或 Issues 来提交代码修改、功能增强或报告问题。

