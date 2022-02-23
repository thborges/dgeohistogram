/*
 * A template to use property like fields and simplify code reading.
 * From https://stackoverflow.com/a/35478633
 * 
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

template <typename T>
class Property {
public:
    virtual ~Property() {}  //C++11: use override and =default;
    virtual T& operator= (const T& f) { return value = f; }
    virtual const T& operator() () const { return value; }
    virtual explicit operator const T& () const { return value; }
    virtual T* operator->() { return &value; }
protected:
    T value;
};

template <typename T>
class ReadOnlyProperty {
public:
    virtual ~ReadOnlyProperty() {}
    virtual operator T const & () const { return value; }
protected:
    T value;
};

/*template <typename T, typename C>
class ReadOnlyPropertyComputed {
public:
    virtual ~ReadOnlyPropertyComputed() {}
    virtual operator T const & () const { return C(); }
protected:
    T value;
};
*/