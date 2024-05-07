# NSGA-III for Chesapeake Bay Watershed BMP Optimization

This repository contains a modified version of the NSGA-III algorithm designed to optimize Best Management Practices (BMP) implementations for the Chesapeake Bay Watershed. The algorithm interfaces with the Chesapeake Assessment Scenario Tool (CAST) to provide environmental planning solutions aimed at reducing nutrient and sediment loads in the watershed.

## Features

- **Integration with CAST**: Utilizes the CoreCAST tool to assess the impact of BMPs on nutrient and sediment loads.
- **Advanced Data Handling**: Leverages AWS S3 for data storage, and utilizes Parquet and Apache Arrow for efficient data processing and storage.

## Dependencies

### System Packages

Ensure the following dependencies are installed:

- AWS SDK for C++
- Parquet (part of Apache Arrow)
- Apache Arrow

You can install these dependencies on most systems using package managers or by compiling from source.

### AWS SDK for C++

```bash
sudo apt-get install libaws-cpp-sdk-s3-dev
```

### Apache Arrow and Parquet

Apache Arrow and Parquet can be installed via Conda or from source. For more information, see the [Apache Arrow installation guide](https://arrow.apache.org/docs/developers/cpp/building.html).

## Installation

To install and compile the modified NSGA-III code, follow these steps:

```bash
git clone https://github.com/gtoscano/nsga3-cbw.git
cd nsga3-cbw/build
cmake ..
make
```

## Usage

After compilation, you can run the program directly from the build directory. Further instructions and parameters can be specified according to your scenario planning needs.

## Contributing

Contributions to improve the algorithm or its integration with CAST are highly appreciated. Please fork this repository and submit pull requests with your suggested changes.

## License

This project is licensed under the Apache License 2.0. For more details, see the [LICENSE](LICENSE) file in the repository.

## Support

For support, queries, or further information, please open an issue in the GitHub repository.

