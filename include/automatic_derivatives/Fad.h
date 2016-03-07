#ifndef Fad_h
#define Fad_h



//############ How to use this class:  ########################
/*

Write this !


// To-do:
//   Move management of value out of UnaryOpFad in Fad_functions.h
//     Note change UnaryOpFad::arg -> UnaryOpFad::x
//   Make power<>(FloatType) forward variables properly
//   Update documentation
// current problems:
//   missing operators: +=, -=, *=, /=  [NEVER FIX!!!]
//   doesn't cancel unnecessary operations (e.g. negation(negation(x)))

*/

//######################### The class definition ##############################
#include <map>
#include <set>
#include <utility>


#define FLOAT_TYPENAME(TypeName) \
typename TypeName, \
  typename = typename std::enable_if<std::is_arithmetic<TypeName>::value>::type

#define FAD_TYPENAME(TypeName) \
typename TypeName, typename = typename std::enable_if<TypeName::isFad>::type


#define FAD_OP_TYPENAME(TypeName) \
typename TypeName, typename = typename std::enable_if<TypeName::isOp>::type


typedef unsigned char Order_T;

template <typename FloatType, Order_T derivOrder>
struct Fad
{
  // typedef a short name for the Fad type
  typedef Fad<FloatType, derivOrder> FadType;
  // typedef a short name for the derivatives' type
  typedef Fad<FloatType, derivOrder-1> DerivType;
  // typedef the float type
  typedef FloatType Float_T;
  
  // default constructor
  Fad();
  // construct and assign value
  Fad(const FloatType & value);
  // construct, assign value and set derivative index
  Fad(const FloatType & value, const size_t & derivIndex);
  // copy constructor
  Fad(const FadType & other);
  // r-value copy constructor
  Fad(FadType && other);
  
  // construct from operation
  template <FAD_OP_TYPENAME(OpType)>
  Fad(const OpType & op) : value (op.value) {
    // copy derivatives
    for(const size_t & index : op.getIndices()) {
      derivatives[index] = op.d(index);
    }
  }
  // construct from operation r-value
  template <FAD_OP_TYPENAME(OpType)>
  Fad(OpType && op) : value (std::move(op.value)) {
    // copy derivatives
    for(const size_t & index : op.getIndices()) {
      derivatives[index] = op.d(index);
    }
  }
  
  // clear the derivative information
  void clear(void);
  
  // set the derivative index for this variable
  void setDerivativeIndex(const size_t & derivIndex);

  // overloading assignment operator
  FadType & operator=(const FadType & other);
  FadType & operator=(FadType && other);
  FadType & operator=(const FloatType & other);
  template <FAD_OP_TYPENAME(OpType)>
  FadType & operator=(const OpType & op) {
    // assign from operation
    
    // copy value
    value = op.value;
    // copy derivatives
    if(derivatives.empty()) {
      // first assignment
      for(const size_t & index : op.getIndices()) {
        derivatives[index] = op.d(index);
      }
    }
    else {
      for(auto & mapObj : derivatives) {
        mapObj.second = op.d(mapObj.first);
      }
    }
  }
  template <FAD_OP_TYPENAME(OpType)>
  FadType & operator=(const OpType && op) {
    // assign from operation r-value
    
    // move value
    value = std::move(op.value);
    // copy derivatives
    if(derivatives.empty()) {
      // first assignment
      for(const size_t & index : op.getIndices()) {
        derivatives[index] = op.d(index);
      }
    }
    else {
      for(auto & mapObj : derivatives) {
        mapObj.second = op.d(mapObj.first);
      }
    }
  }
  
  // check if an index has a nonzero derivative
  bool hasD(const size_t & index) const;
  
  // return the derivative of the Fad without checking indices first
  const DerivType & D(const size_t & index) const;

  // return the derivative of the Fad as another Fad object
  DerivType d(const size_t & index) const;
  
  // return derivOrder
  static constexpr Order_T getDerivOrder(void) { return derivOrder; }
  
  // return the indices to non-zero derivatives
  std::set<size_t> getIndices(void) const;
  
  // expose that this is a Fad-type
  static const bool isFad = true;

  // non-differentiated value
  FloatType value;
  
  // where derivatives are stored:
  std::map<size_t,DerivType> derivatives;
};



template <typename FloatType>
struct Fad<FloatType, 1>
{
  // typedef a short name for the Fad type
  typedef Fad<FloatType, 1> FadType;
  // typedef a short name for the derivatives' type
  typedef FloatType DerivType;
  // typedef the float type
  typedef FloatType Float_T;
    
   // default constructor
  Fad();
  // construct and assign value
  Fad(const FloatType & value);
  // construct, assign value and set derivative index
  Fad(const FloatType & value, const size_t & derivIndex);
  // copy constructor
  Fad(const FadType & other);
  // r-value copy constructor
  Fad(FadType && other);
  
  // construct from operation
  template <FAD_OP_TYPENAME(OpType)>
  Fad(const OpType & op) : value (op.value) {
    // copy derivatives
    for(const size_t & index : op.getIndices()) {
      derivatives[index] = op.D(index);
    }
  }
  // construct from operation r-value
  template <FAD_OP_TYPENAME(OpType)>
  Fad(OpType && op) : value (std::move(op.value)) {
    // copy derivatives
    for(const size_t & index : op.getIndices()) {
      derivatives[index] = op.D(index);
    }
  }
  
  // clear the derivative information
  void clear(void);
  
  // set the derivative index for this variable
  void setDerivativeIndex(const size_t & derivIndex);

  // overloading assignment operator
  FadType & operator=(const FadType & other);
  FadType & operator=(FadType && other);
  FadType & operator=(const FloatType & other);
  template <FAD_OP_TYPENAME(OpType)>
  FadType & operator=(const OpType & op) {
    // assign from operation
    
    // copy value
    value = op.getValue();
    // copy derivatives
    if(derivatives.empty()) {
      // first assignment
      for(const size_t & index : op.getIndices()) {
        derivatives[index] = op.D(index);
      }
    }
    else {
      for(auto & mapObj : derivatives) {
        mapObj.second = op.D(mapObj.first);
      }
    }
    return *this;
  }
  template <FAD_OP_TYPENAME(OpType)>
  FadType & operator=(const OpType && op) {
    // assign from operation r-value
    
    // move value
    value = std::move(op.value);
    // copy derivatives
    if(derivatives.empty()) {
      // first assignment
      for(const size_t & index : op.getIndices()) {
        derivatives[index] = op.D(index);
      }
    }
    else {
      for(auto & mapObj : derivatives) {
        mapObj.second = op.D(mapObj.first);
      }
    }
    return *this;
  }
  
  // check if an index has a nonzero derivative
  bool hasD(const size_t & index) const;
  
  // return the derivative of the Fad without checking indices first
  const DerivType & D(const size_t & index) const;

  // return the derivative of the Fad as another Fad object
  DerivType d(const size_t & index) const;
  
  // return derivOrder
  static constexpr Order_T getDerivOrder(void) { return 1; }
  
  // return the indices to non-zero derivatives
  std::set<size_t> getIndices(void) const;
  
  // expose that this is a Fad-type
  static const bool isFad = true;

  // non-differentiated value
  FloatType value;
    
  // where derivatives are stored:
  std::map<size_t,FloatType> derivatives;
};


//################## Assigning and constructing Fad objects ###################
#include "Fad_construct.h"

//################## Arithmetic operators and helper classes ##################
#include "Fad_operators.h"

//################### special functions (sin/cos/etc) #########################
#include "Fad_functions.h"

#endif
