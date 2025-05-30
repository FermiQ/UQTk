# Array1D.h

## Overview

The `Array1D.h` file defines a C++ template class `Array1D<T>` designed to store and manage one-dimensional arrays of any data type `T`. It provides a dynamic array implementation with functionalities similar to `std::vector` but includes additional features such as Fortran-like parenthesis `()` for element access, methods for direct pointer access (useful for C/Fortran interoperability), and binary I/O operations. The file also includes explicit template specializations for `Array1D<int>` and `Array1D<double>` which mirror the generic template's functionality but are optimized for these common types.

## Key Components

*   **`Array1D<T>` (Generic Template Class)**:
    *   Manages a 1D array of elements of type `T`.
    *   Provides constructors for creating empty arrays, arrays of a specific size, or arrays initialized with a specific value.
    *   Offers methods for resizing, clearing, accessing elements (via `()` and `[]` operators), getting the size/length, and obtaining raw data pointers.
    *   Includes methods for inserting/erasing elements, pushing elements to the back, and setting all values to a specific value.
    *   Supports binary serialization (dumping to/reading from files).
    *   Contains methods specifically for Python interfacing (`shape`, `assign`, `setArray`, `flatten`, `type`, `DumpBinary4py`, `ReadBinary4py`).
    *   Member variables `xsize_` (size) and `data_` (std::vector storing elements) are public for easier Python wrapping.

*   **`Array1D<int>` (Template Specialization)**:
    *   A specialized version of `Array1D` for `int` data type.
    *   It largely replicates the public interface and functionality of the generic `Array1D<T>` template.
    *   Includes specific methods `setnpintArray` and `getnpintArray` for interaction with NumPy-like integer arrays (via pointers/vectors).

*   **`Array1D<double>` (Template Specialization)**:
    *   A specialized version of `Array1D` for `double` data type.
    *   It also largely replicates the public interface and functionality of the generic `Array1D<T>` template.
    *   Includes specific methods `setnpdblArray` and `getnpdblArray` for interaction with NumPy-like double arrays (via pointers/vectors).

## Important Variables/Constants

For `Array1D<T>` (and its specializations):
*   **`int xsize_`**: Public member variable storing the number of elements in the array.
*   **`std::vector<T> data_`**: Public member variable (a `std::vector`) that holds the actual array elements. Making this public facilitates easier integration with Python.

## Usage Examples

```cpp
#include "Array1D.h" // Assuming this is in the include path
#include <iostream>
#include <string>
#include <vector>

int main() {
    // Using the generic template Array1D<T> with doubles
    Array1D<double> dblArray(5, 0.0); // Create an array of 5 doubles, initialized to 0.0
    dblArray(0) = 1.1;
    dblArray[1] = 2.2; // Can also use [] operator
    dblArray.PushBack(3.3); // Add an element to the end

    std::cout << "Double Array (size " << dblArray.XSize() << "):" << std::endl;
    for (int i = 0; i < dblArray.Length(); ++i) {
        std::cout << "dblArray(" << i << ") = " << dblArray(i) << std::endl;
    }

    // Using the Array1D<int> specialization
    Array1D<int> intArray; // Default constructor
    intArray.Resize(3, 10); // Resize to 3 elements, all initialized to 10
    intArray.insert(20, 1); // Insert 20 at index 1

    std::cout << "\nInt Array (size " << intArray.XSize() << "):" << std::endl;
    std::vector<int> intVec = intArray.flatten(); // Get data as std::vector
    for (size_t i = 0; i < intVec.size(); ++i) {
        std::cout << "intArray[" << i << "] = " << intVec[i] << std::endl;
    }

    // Example with std::string
    Array1D<std::string> strArray(2, "hello");
    strArray.PushBack("world");
    std::cout << "\nString Array (size " << strArray.XSize() << "):" << std::endl;
    for (int i = 0; i < strArray.XSize(); ++i) {
        std::cout << strArray(i) << " ";
    }
    std::cout << std::endl;
    
    // Getting array pointer (use with caution)
    if (dblArray.XSize() > 0) {
        double* pDbl = dblArray.GetArrayPointer();
        pDbl[0] = 100.1; // Modify directly
        std::cout << "\nFirst element of dblArray after pointer modification: " << dblArray(0) << std::endl;
    }

    return 0;
}

```

## Dependencies and Interactions

*   **Standard Library Headers**:
    *   `<string>`, `<string.h>`: For string manipulation (though `string.h` is C-style, `<string>` is C++). `type()` method in generic template returns "string".
    *   `<iostream>`: For I/O operations (e.g., `DumpBinary4py`, `ReadBinary4py` use `ofstream`/`ifstream`).
    *   `<vector>`: The internal data storage (`data_`) is an `std::vector<T>`.
    *   `<fstream>`: Used for file I/O operations (`DumpBinary4py`, `ReadBinary4py`).
    *   `<iterator>`, `<algorithm>`: Potentially used by `std::vector` operations or for methods like `copy` in specializations.
    *   `<typeinfo>`: Potentially for type introspection, though not explicitly used in the public interface shown.
*   **`error_handlers.h`**: Included for UQTk's error handling mechanisms (e.g., `Tantrum` exception, though commented out in `insert`/`erase` methods in the provided header snippet).
*   **C/Fortran Interoperability**: The `GetArrayPointer()` and `GetConstArrayPointer()` methods are designed to facilitate passing the array data to C or Fortran routines that expect raw pointers.
*   **Python Interfacing**: Several methods are explicitly commented as being for Python interfacing (e.g., `shape()`, `assign()`, `flatten()`, `type()`, `DumpBinary4py()`, `ReadBinary4py()`, `setnpintArray()`, `getnpintArray()`, `setnpdblArray()`, `getnpdblArray()`). The public nature of `xsize_` and `data_` also supports this.
```
