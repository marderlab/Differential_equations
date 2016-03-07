#ifndef Fad_h
#error Do not include this file directly, instead include Fad.h
#endif
#ifndef Fad_operators_h
#define Fad_operators_h

// This file defines basic relational and arithmetic operators:
//  < > <= >= == != + -(negation) -(subtraction) * /
// To do arithmetic, helper classes are required (behind the scenes):
// WrapFad, AddFad, NegFad, SubFad, MulFad, PowFad, DivFad


//////////////////////////// Helper functions /////////////////////////////////

template <class T>
std::set<T> setUnion(std::set<T> a, std::set<T> b)
{
  // return the union of sets a and b
  
  for(const T & b_n : b) {
    // insert the elements of b into set a
    a.insert(std::move(b_n));
  }
  // return the enlarged set a
  return a;
}



////////////////////////// Relational operators  //////////////////////////////

#define DEFINE_FAD_RELATIONAL_OPERATOR(op) \
template<typename FadType, \
         typename = typename std::enable_if<FadType::isFad>::type> \
inline bool operator op(typename FadType::Float_T x, \
                          const FadType & y) \
{ \
  return x op y.value; \
} \
template<typename FadType, \
         typename = typename std::enable_if<FadType::isFad>::type> \
inline bool operator op(const FadType & x, \
                      typename FadType::Float_T y) \
{ \
  return x.value op y; \
} \
template<typename XFadType, \
         typename YFadType, \
         typename = typename std::enable_if<YFadType::isFad && \
                                            XFadType::isFad>::type> \
inline bool operator op(const XFadType & x, \
                      const YFadType & y) \
{ \
  return x.value op y.value; \
}

DEFINE_FAD_RELATIONAL_OPERATOR(<)
DEFINE_FAD_RELATIONAL_OPERATOR(>)
DEFINE_FAD_RELATIONAL_OPERATOR(<=)
DEFINE_FAD_RELATIONAL_OPERATOR(>=)
DEFINE_FAD_RELATIONAL_OPERATOR(==)
DEFINE_FAD_RELATIONAL_OPERATOR(!=)

#undef DEFINE_FAD_RELATIONAL_OPERATOR



//////////////////// OpBase: basis struct for all operators ///////////////////

template <typename FloatType, Order_T derivOrder>
struct OpBase
{
  // typedef the float type
  typedef FloatType Float_T;
  
  // return derivOrder
  static constexpr Order_T getDerivOrder(void) { return derivOrder; }

  // expose that this is an operator
  static const bool isOp = true;
  // expose that this is a Fad-type
  static const bool isFad = true;
};



//////////// RefArgType: struct for managing if an arg is a reference /////////

template <typename ArgType, bool refArg>
struct RefArgType
{
  typedef const ArgType type;
};
template <typename ArgType>
struct RefArgType<ArgType, true>
{
  typedef const ArgType & type;
};



//////////////////// UnaryOpFad: struct for unary operators ///////////////////

template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct UnaryOpFad : public OpBase<FloatType, derivOrder>
{
  // declare x
  typename RefArgType<ArgType, refArg>::type x;

  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  UnaryOpFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      x (std::move(arg)) {}
  template<bool U = refArg>
  UnaryOpFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      x (arg) {}
  
  // check if an index has a nonzero derivative
  bool hasD(const size_t & index) const { return x.hasD(index); }
  
  // return the indices to non-zero derivatives
  std::set<size_t> getIndices(void) const { return x.getIndices(); }
};



/////////////////////////////// wrap function /////////////////////////////////

// declare WrapFad struct (define it later)
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct WrapFad;

template <Order_T derivOrder, FAD_TYPENAME(ArgType)>
inline WrapFad<typename ArgType::Float_T, ArgType, derivOrder, false>
wrap(ArgType && x)
{
  return WrapFad<typename ArgType::Float_T, ArgType, derivOrder, false>
    (std::move(x));
}
template <Order_T derivOrder, FAD_TYPENAME(ArgType)>
inline WrapFad<typename ArgType::Float_T, ArgType, derivOrder, true>
wrap(const ArgType & x)
{
  return WrapFad<typename ArgType::Float_T, ArgType, derivOrder, true> (x);
}



/////////////////////////// Helper function: getD /////////////////////////////
// auto getD<derivOrder>(arg, index) gets the derivative of arg packaged in
// a Fad-type of order derivOrder

template <Order_T derivOrder, typename ArgType,
          typename = typename std::enable_if<derivOrder != 1 &&
                                  ArgType::getDerivOrder() == derivOrder>::type>
auto
getD(const ArgType & arg, const size_t & index)
  -> decltype(arg.d(index))
{
  return arg.d(index);
}
template <Order_T derivOrder, typename ArgType,
          typename = typename std::enable_if<derivOrder != 1 &&
                                  ArgType::getDerivOrder() != derivOrder>::type>
auto
getD(const ArgType & arg, const size_t & index)
  -> decltype(wrap<derivOrder>(arg.d(index)))
{
  return wrap<derivOrder>(arg.d(index));
}
template <Order_T derivOrder, typename ArgType,
          typename = typename std::enable_if<derivOrder == 1 &&
                                           ArgType::getDerivOrder() != 1>::type>
auto
getD(const ArgType & arg, const size_t & index)
  -> decltype(arg.d(index).value)
{
  return arg.d(index).value;
}
template <Order_T derivOrder, typename ArgType,
          typename = typename std::enable_if<derivOrder == 1 &&
                                           ArgType::getDerivOrder() == 1>::type>
auto getD(const ArgType & arg, const size_t & index)
  -> decltype(arg.D(index))
{
  return arg.D(index);
}



//////////////////// WrapFad - class for adjusting derivOrder /////////////////

template <typename FloatType, typename ArgType, Order_T derivOrder, bool refArg>
struct WrapFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;

  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  WrapFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      value (x.value) {}
  template<bool U = refArg>
  WrapFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
     value (x.value) {}
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const -> decltype(getD<derivOrder>(x, index)){
    return getD<derivOrder>(x, index);
  }
  
  typedef decltype(((WrapFad*)NULL)->d(0)) DerivType;
  const FloatType & value;
};

template <typename FloatType, typename ArgType, bool refArg>
struct WrapFad<FloatType, ArgType, 1, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef a short name for the derivatives' type
  typedef FloatType DerivType;
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  WrapFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      value (x.value) {}
  template<bool U = refArg>
  WrapFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
     value (x.value) {}
  
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return x.hasD(index) ? getD<1>(x, index) : (FloatType)0.0;
  }
  auto D(const size_t & index) const -> decltype(getD<1>(x, index))
  {
    return getD<1>(x, index);
  }
  
  const FloatType & value;
};



//////////////////////////// NegFad - Negation class ///////////////////////////

// -Fad with derivOrder > 1
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct NegFad : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  NegFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)), value(-x.value) {}
  template<bool U = refArg>
  NegFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op), value (-x.value) {}
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(-getD<derivOrder>(x, index)){
    return -getD<derivOrder>(x, index);
  }
  
  typedef decltype(((NegFad*)NULL)->d(0)) DerivType;
  
  const FloatType value;
};


// -Fad  with derivOrder == 1, Fad order > 1
template <typename FloatType, typename ArgType, bool refArg>
struct NegFad<FloatType, ArgType, 1, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  NegFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)), value (-x.value) {}
  template<bool U = refArg>
  NegFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op), value(-x.value) {}

  // return the derivative as a FloatType
  FloatType d(const size_t & index) const {
    return x.hasD(index) ? -getD<1>(x, index) : (FloatType)0.0;
  }
  
  // return the derivative as a FloatType
  FloatType D(const size_t & index) const {
    return -getD<1>(x, index);
  }

  // typedef a short name for the derivatives' type
  typedef FloatType DerivType;
  
  const FloatType value;
};



//////////////////////////// Negation operator ////////////////////////////////

template <FAD_TYPENAME(ArgType)>
inline NegFad<typename ArgType::Float_T, ArgType, ArgType::getDerivOrder(),
              false>
operator-(ArgType && x)
{
  return NegFad<typename ArgType::Float_T, ArgType, ArgType::getDerivOrder(),
                false> (std::move(x));
}
template <FAD_TYPENAME(ArgType)>
inline NegFad<typename ArgType::Float_T, ArgType, ArgType::getDerivOrder(),
              true>
operator-(const ArgType & x)
{
  return NegFad<typename ArgType::Float_T, ArgType, ArgType::getDerivOrder(),
                true> (x);
}



/////////////////// BinaryOpFad: struct for binary operators //////////////////


template <typename FloatType, typename LeftType, typename RightType,
          Order_T derivOrder, bool refLeft, bool refRight>
struct BinaryOpFad : public OpBase<FloatType, derivOrder>
{
  // define constructor, with and without pass-by-reference
  template<bool L = refLeft, bool R = refRight>
  BinaryOpFad(LeftType && op1, RightType && op2,
              typename std::enable_if<!L&!R>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      left (std::move(op1)),
      right (std::move(op2)) {}
  template<bool L = refLeft, bool R = refRight>
  BinaryOpFad(LeftType && op1, const RightType & op2,
              typename std::enable_if<!L&R>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      left (std::move(op1)),
      right (op2) {}
  template<bool L = refLeft, bool R = refRight>
  BinaryOpFad(const LeftType & op1, RightType && op2,
              typename std::enable_if<L&!R>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      left (op1),
      right (std::move(op2)) {}
  template<bool L = refLeft, bool R = refRight>
  BinaryOpFad(const LeftType & op1, const RightType & op2,
              typename std::enable_if<L&R>::type* = 0)
    : OpBase<FloatType, derivOrder>(),
      left (op1),
      right (op2) {}

  // check if an index has a nonzero derivative
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<!LFloat&!RFloat, bool>::type
  hasD(const size_t & index) const {
    return left.hasD(index) || right.hasD(index);
  }
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<!LFloat&RFloat, bool>::type
  hasD(const size_t & index) const {
    return left.hasD(index);
  }
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<LFloat&!RFloat, bool>::type
  hasD(const size_t & index) const {
    return right.hasD(index);
  }
  
  // return the indices to non-zero derivatives
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<!LFloat&!RFloat, std::set<size_t>>::type
  getIndices(void) const {
    return setUnion(left.getIndices(), right.getIndices());
  }
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<!LFloat&RFloat, std::set<size_t>>::type
  getIndices(void) const {
    return left.getIndices();
  }
  template<bool LFloat = std::is_same<LeftType, FloatType>::value,
           bool RFloat = std::is_same<RightType, FloatType>::value>
  typename std::enable_if<LFloat&!RFloat, std::set<size_t>>::type
  getIndices(void) const {
    return right.getIndices();
  }

  typename RefArgType<LeftType, refLeft>::type left;
  typename RefArgType<RightType, refRight>::type right;
};


///////// Macros to define constructor of BinaryOpFad-based operators /////////

#define BINARY_CONSTRUCTOR(className, LType, RType, dOrder, initValue) \
  typedef BinaryOpFad<FloatType, LType, RType, dOrder, refLeft, refRight> \
    OpType; \
  using OpType::left;  using OpType::right; \
  const FloatType value; \
  \
  template<bool L = refLeft, bool R = refRight> \
  className (LType && op1, RType && op2, \
         typename std::enable_if<!L&!R>::type* = 0) \
    : OpType(std::move(op1), std::move(op2)), value (initValue) {} \
  template<bool L = refLeft, bool R = refRight> \
  className (LType && op1, const RType & op2, \
         typename std::enable_if<!L&R>::type* = 0) \
    : OpType(std::move(op1), op2), value (initValue) {} \
  template<bool L = refLeft, bool R = refRight> \
  className (const LType & op1, RType && op2, \
        typename std::enable_if<L&!R>::type* = 0) \
    : OpType(op1, std::move(op2)), value (initValue) {} \
  template<bool L = refLeft, bool R = refRight> \
  className (const LType & op1, const RType & op2, \
         typename std::enable_if<L&R>::type* = 0) \
    : OpType(op1, op2), value (initValue) {}


/////////////////// Macro for creating arithmetic operators ///////////////////

// when performing Fad op Fad, return an object with order equal to the minimum
// of the two operands
template <typename LeftType, typename RightType>
constexpr Order_T minOrder(void)
{
  return LeftType::getDerivOrder() < RightType::getDerivOrder() ?
            LeftType::getDerivOrder() : RightType::getDerivOrder();
}

// This macro defines operator functions for the following pairs:
//   Fad op Fad  (by const reference and r-value)
//   Fad op Fad::Float_T  (by const reference and r-value)
//   Fad::Float_T op Fad  (by const reference and r-value)
//   Fad op FloatType  (Fad by const reference and r-value, FloatType by value)
//   FloatType op Fad  (Fad by const reference and r-value, FloatType by value)
#define DEFINE_FAD_ARITHMETIC_OPERATOR(opName, op) \
template <typename LeftType, typename RightType, \
          typename = typename std::enable_if<LeftType::isFad && \
                                             RightType::isFad>::type> \
inline auto \
operator op (LeftType && x, RightType && y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
            minOrder<LeftType, RightType>(), \
            false, false> \
{ \
  static_assert(std::is_same<typename LeftType::Float_T, \
                             typename RightType::Float_T>::value, \
                "x and y have incompatible float types"); \
  return opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
                minOrder<LeftType, RightType>(), \
                false, false> \
         (std::move(x), std::move(y)); \
} \
template <typename LeftType, typename RightType, \
          typename = typename std::enable_if<LeftType::isFad && \
                                             RightType::isFad>::type> \
inline auto \
operator op (LeftType && x, const RightType & y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
            minOrder<LeftType, RightType>(), \
            false, true> \
{ \
  static_assert(std::is_same<typename LeftType::Float_T, \
                             typename RightType::Float_T>::value, \
                "x and y have incompatible float types"); \
  return opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
                minOrder<LeftType, RightType>(), \
                false, true> \
         (std::move(x), y); \
} \
template <typename LeftType, typename RightType, \
          typename = typename std::enable_if<LeftType::isFad && \
                                             RightType::isFad>::type> \
inline auto \
operator op (const LeftType & x, RightType && y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
            minOrder<LeftType, RightType>(), \
            true, false> \
{ \
  static_assert(std::is_same<typename LeftType::Float_T, \
                             typename RightType::Float_T>::value, \
                "x and y have incompatible float types"); \
  return opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
                minOrder<LeftType, RightType>(), \
                true, false> \
         (x, std::move(y)); \
} \
template <typename LeftType, typename RightType, \
          typename = typename std::enable_if<LeftType::isFad && \
                                             RightType::isFad>::type> \
inline auto \
operator op (const LeftType & x, const RightType & y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
            minOrder<LeftType, RightType>(), \
            true, true> \
{ \
  static_assert(std::is_same<typename LeftType::Float_T, \
                             typename RightType::Float_T>::value, \
                "x and y have incompatible float types"); \
  return opName##Fad<typename LeftType::Float_T, LeftType, RightType, \
                minOrder<LeftType, RightType>(), \
                true, true> \
         (x, y); \
} \
\
\
template <typename LeftType, \
          typename = typename std::enable_if<LeftType::isFad>::type>\
inline auto \
operator op (LeftType && x, typename LeftType::Float_T && y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), false, false> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                false, false> \
           ( std::move(x), std::move(y) ); \
} \
template <typename LeftType, \
          typename = typename std::enable_if<LeftType::isFad>::type>\
inline auto \
operator op (const LeftType & x, typename LeftType::Float_T && y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), true, false> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                true, false> \
           ( x, std::move(y) ); \
} \
template <typename LeftType, \
          typename = typename std::enable_if<LeftType::isFad>::type>\
inline auto \
operator op (LeftType && x, const typename LeftType::Float_T & y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), false, true> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                false, true> \
           ( std::move(x), y ); \
} \
template <typename LeftType, \
          typename = typename std::enable_if<LeftType::isFad>::type>\
inline auto \
operator op (const LeftType & x, const typename LeftType::Float_T & y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), true, true> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                true, true> \
           ( x, y ); \
} \
\
template <typename LeftType, typename FloatType, \
          typename = typename std::enable_if<LeftType::isFad && \
            std::is_arithmetic<FloatType>::value && \
            !std::is_same<FloatType, typename LeftType::Float_T>::value>::type>\
inline auto \
operator op (LeftType && x, FloatType y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), false, false> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                false, false> \
           ( std::move(x), (typename LeftType::Float_T)y ); \
} \
template <typename LeftType, typename FloatType, \
          typename = typename std::enable_if<LeftType::isFad && \
            std::is_arithmetic<FloatType>::value && \
            !std::is_same<FloatType, typename LeftType::Float_T>::value>::type>\
inline auto \
operator op (const LeftType & x, FloatType y) \
  -> opName##Fad<typename LeftType::Float_T, LeftType, typename LeftType::Float_T, \
            LeftType::getDerivOrder(), true, false> \
{ \
  return opName##Fad<typename LeftType::Float_T, LeftType, \
                typename LeftType::Float_T, LeftType::getDerivOrder(), \
                true, false> \
           ( x, (typename LeftType::Float_T)(y) ); \
} \
\
\
template <typename RightType, \
        typename = typename std::enable_if<RightType::isFad>::type>\
inline auto \
operator op (typename RightType::Float_T && x, RightType && y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), false, false> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
                RightType, RightType::getDerivOrder(), false, false> \
           ( std::move(x), std::move(y) ); \
} \
template <typename RightType, \
        typename = typename std::enable_if<RightType::isFad>::type>\
inline auto \
operator op (typename RightType::Float_T && x, const RightType & y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), false, true> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
                RightType, RightType::getDerivOrder(), false, true> \
           ( std::move(x), y ); \
} \
template <typename RightType, \
        typename = typename std::enable_if<RightType::isFad>::type>\
inline auto \
operator op (const typename RightType::Float_T & x, RightType && y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), true, false> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
                RightType, RightType::getDerivOrder(), true, false> \
           ( x, std::move(y) ); \
} \
template <typename RightType, \
         typename = typename std::enable_if<RightType::isFad>::type> \
inline auto \
operator op (const typename RightType::Float_T & x, const RightType & y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), true, true> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
                RightType, RightType::getDerivOrder(), true, true> \
           ( x, y ); \
} \
\
template <typename FloatType, typename RightType, \
          typename = typename std::enable_if<RightType::isFad && \
            std::is_arithmetic<FloatType>::value && \
            !std::is_same<FloatType, typename RightType::Float_T>::value>::type> \
inline auto \
operator op (FloatType x, RightType && y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), false, false> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
                RightType, RightType::getDerivOrder(), false, false> \
           ( (typename RightType::Float_T)x, std::move(y) ); \
} \
template <typename FloatType, typename RightType, \
          typename = typename std::enable_if<RightType::isFad && \
            std::is_arithmetic<FloatType>::value && \
            !std::is_same<FloatType, typename RightType::Float_T>::value>::type> \
inline auto \
operator op (FloatType x, const RightType & y) \
  -> opName##Fad<typename RightType::Float_T, typename RightType::Float_T, \
            RightType, RightType::getDerivOrder(), false, true> \
{ \
  return opName##Fad<typename RightType::Float_T, typename RightType::Float_T,\
                RightType, RightType::getDerivOrder(), false, true> \
           ( (typename RightType::Float_T)x, y ); \
}



//////////////////////////// AddFad - Addition class //////////////////////////

// Fad + Fad with derivOrder > 1
template <typename FloatType, typename LeftType, typename RightType,
           Order_T derivOrder, bool refLeft=false, bool refRight=false>
struct AddFad
  : public BinaryOpFad<FloatType, LeftType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, LeftType, RightType, derivOrder,
                     left.value + right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(left, index) + getD<derivOrder>(right, index))
  {
    return getD<derivOrder>(left, index) + getD<derivOrder>(right, index);
  }
  
  typedef decltype(((AddFad*)NULL)->d(0)) DerivType;
};

// Fad + Float with derivOrder > 1
template <typename FloatType, typename LeftType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct AddFad<FloatType, LeftType, FloatType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, LeftType, FloatType, derivOrder,
                     left.value + right)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const -> decltype(getD<derivOrder>(left, index))
  {
    return getD<derivOrder>(left, index);
  }
  
  typedef decltype(((AddFad*)NULL)->d(0)) DerivType;
};

// Float + Fad with derivOrder > 1
template <typename FloatType, typename RightType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct AddFad<FloatType, FloatType, RightType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, FloatType, RightType, derivOrder,
                     left + right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const -> decltype(getD<derivOrder>(right, index))
  {
    return getD<derivOrder>(right, index);
  }
  
  typedef decltype(((AddFad*)NULL)->d(0)) DerivType;
};

// Fad + Fad with derivOrder == 1
template <typename FloatType, typename LeftType, typename RightType,
          bool refLeft, bool refRight>
struct AddFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, LeftType, RightType, 1, left.value + right.value)

  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? right.hasD(index) ? getD<1>(left, index)
                                                  + getD<1>(right, index)
                                                : getD<1>(left, index)
                            : right.hasD(index) ? getD<1>(right, index)
                                                : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return left.hasD(index) ? right.hasD(index) ? getD<1>(left, index)
                                                  + getD<1>(right, index)
                                                : getD<1>(left, index)
                            : getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};

// Fad + Float with derivOrder == 1
template <typename FloatType, typename LeftType, bool refLeft, bool refRight>
struct AddFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, LeftType, FloatType, 1, left.value + right)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? getD<1>(left, index) : (FloatType)0.0;
  }
  auto D(const size_t & index) const -> decltype(getD<1>(left, index))
  {
    return getD<1>(left, index);
  }
  
  typedef FloatType DerivType;
};

// Float + Fad with derivOrder == 1
template <typename FloatType, typename RightType, bool refLeft, bool refRight>
struct AddFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(AddFad, FloatType, RightType, 1, left + right.value)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return right.hasD(index) ? getD<1>(right, index) : (FloatType)0.0;
  }
  auto D(const size_t & index) const -> decltype(getD<1>(right, index))
  {
    return getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};


DEFINE_FAD_ARITHMETIC_OPERATOR(Add, +)



////////////////////////// SubFad - Subtraction class //////////////////////////

// Fad - Fad with derivOrder > 1
template <typename FloatType, typename LeftType, typename RightType,
           Order_T derivOrder, bool refLeft=false, bool refRight=false>
struct SubFad
  : public BinaryOpFad<FloatType, LeftType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, LeftType, RightType, derivOrder,
                     left.value - right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(left, index) - getD<derivOrder>(right, index))
  {
    return getD<derivOrder>(left, index) - getD<derivOrder>(right, index);
  }
  
  typedef decltype(((SubFad*)NULL)->d(0)) DerivType;
};

// Fad - Float with derivOrder > 1
template <typename FloatType, typename LeftType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct SubFad<FloatType, LeftType, FloatType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, LeftType, FloatType, derivOrder,
                     left.value - right)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const -> decltype(getD<derivOrder>(left, index))
  {
    return getD<derivOrder>(left, index);
  }
  
  typedef decltype(((SubFad*)NULL)->d(0)) DerivType;
};

// Float - Fad with derivOrder > 1
template <typename FloatType, typename RightType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct SubFad<FloatType, FloatType, RightType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, FloatType, RightType, derivOrder,
                     left - right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const -> decltype(-getD<derivOrder>(right, index))
  {
    return -getD<derivOrder>(right, index);
  }
  
  typedef decltype(((SubFad*)NULL)->d(0)) DerivType;
};

// Fad - Fad with derivOrder == 1
template <typename FloatType, typename LeftType, typename RightType,
          bool refLeft, bool refRight>
struct SubFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, LeftType, RightType, 1, left.value - right.value)

  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? right.hasD(index) ? getD<1>(left, index)
                                                  - getD<1>(right, index)
                                                : getD<1>(left, index)
                            : right.hasD(index) ? -getD<1>(right, index)
                                                : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return left.hasD(index) ? right.hasD(index) ? getD<1>(left, index)
                                                  - getD<1>(right, index)
                                                : getD<1>(left, index)
                            : -getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};

// Fad - Float with derivOrder == 1
template <typename FloatType, typename LeftType, bool refLeft, bool refRight>
struct SubFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, LeftType, FloatType, 1, left.value - right)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? getD<1>(left, index) : (FloatType)0.0;
  }
  auto D(const size_t & index) const -> decltype(getD<1>(left, index))
  {
    return getD<1>(left, index);
  }
  
  typedef FloatType DerivType;
};

// Float - Fad with derivOrder == 1
template <typename FloatType, typename RightType, bool refLeft, bool refRight>
struct SubFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(SubFad, FloatType, RightType, 1, left - right.value)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return right.hasD(index) ? -getD<1>(right, index) : (FloatType)0.0;
  }
  auto D(const size_t & index) const -> decltype(-getD<1>(right, index))
  {
    return -getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};


DEFINE_FAD_ARITHMETIC_OPERATOR(Sub, -)



///////////////////////// MulFad - Multiplication class ///////////////////////

// Fad * Fad with derivOrder > 1
template <typename FloatType, typename LeftType, typename RightType,
           Order_T derivOrder, bool refLeft=false, bool refRight=false>
struct MulFad
  : public BinaryOpFad<FloatType, LeftType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, LeftType, RightType, derivOrder,
                     left.value * right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(left, index) * right
                + left * getD<derivOrder>(right, index))
  {
    return getD<derivOrder>(left, index) * right
         + left * getD<derivOrder>(right, index);
  }
  
  typedef decltype(((MulFad*)NULL)->d(0)) DerivType;
};

// Fad * Float with derivOrder > 1
template <typename FloatType, typename LeftType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct MulFad<FloatType, LeftType, FloatType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, LeftType, FloatType, derivOrder,
                     left.value * right)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(left, index) * right)
  {
    return getD<derivOrder>(left, index) * right;
  }
  
  typedef decltype(((MulFad*)NULL)->d(0)) DerivType;
};

// Float * Fad with derivOrder > 1
template <typename FloatType, typename RightType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct MulFad<FloatType, FloatType, RightType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, FloatType, RightType, derivOrder,
                     left * right.value)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(left * getD<derivOrder>(right, index))
  {
    return left * getD<derivOrder>(right, index);
  }
  
  typedef decltype(((MulFad*)NULL)->d(0)) DerivType;
};

// Fad * Fad with derivOrder == 1
template <typename FloatType, typename LeftType, typename RightType,
          bool refLeft, bool refRight>
struct MulFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, LeftType, RightType, 1, left.value * right.value)

  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index)
      ? right.hasD(index) ? getD<1>(left, index) * right.value
                            + left.value * getD<1>(right, index)
                          : getD<1>(left, index) * right.value
      : right.hasD(index) ? left.value * getD<1>(right, index)
                          : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return left.hasD(index)
      ? right.hasD(index) ? getD<1>(left, index) * right.value
                            + left.value * getD<1>(right, index)
                          : getD<1>(left, index) * right.value
      : left.value * getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};

// Fad * Float with derivOrder == 1
template <typename FloatType, typename LeftType, bool refLeft, bool refRight>
struct MulFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, LeftType, FloatType, 1, left.value * right)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? getD<1>(left, index) * right : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return getD<1>(left, index) * right;
  }
  
  typedef FloatType DerivType;
};

// Float * Fad with derivOrder == 1
template <typename FloatType, typename RightType, bool refLeft, bool refRight>
struct MulFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(MulFad, FloatType, RightType, 1, left * right.value)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return right.hasD(index) ? left * getD<1>(right, index) : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return left * getD<1>(right, index);
  }
  
  typedef FloatType DerivType;
};


DEFINE_FAD_ARITHMETIC_OPERATOR(Mul, *)



/////////////////// power: compile-time integer exponentiation ////////////////

// power must be known at compile time, e.g. to cube
// use as xCubed = power<3>(x), where x can be any type

template <size_t p, typename T, bool parity>
struct PowClass;

template <size_t p, typename T>
struct PowClass<p, T, true>
{
  auto power(const T & x) const
    -> decltype(x * PowClass<p/2, T, (p/2) % 2>().power(x)
                  * PowClass<p/2, T, (p/2) % 2>().power(x)) {
    auto y = PowClass<p/2, T, (p/2) % 2>().power(x);
    return x * y * y;
  }
};
template <size_t p, typename T>
struct PowClass<p, T, false>
{
  auto power(const T & x) const
    -> decltype(PowClass<p/2, T, (p/2) % 2>().power(x)
              * PowClass<p/2, T, (p/2) % 2>().power(x)) {
    auto y = PowClass<p/2, T, (p/2) % 2>().power(x);
    return y * y;
  }
};
template <class T>
struct PowClass<1, T, true>
{
  T power(const T & x) const {
    return x;
  }
};
template <class T>
struct PowClass<0, T, false>
{
  T power(void) const {
    return 1;
  }
};

template <size_t p, FLOAT_TYPENAME(T)>
inline T power(const T & x)
{
  return PowClass<p, T, p % 2>().power(x);
}



//////////////////////// PowFad - Exponentiation class ////////////////////////

// -Fad with derivOrder > 1, p > 1
template <typename FloatType, typename ArgType, Order_T derivOrder,
          size_t p, bool refArg=false>
struct PowFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the power less one
  typedef PowFad<FloatType, ArgType, derivOrder-1, p-1, true> PowM1Type;
  // typedef the derivative formula
  typedef MulFad<FloatType, FloatType, PowM1Type, derivOrder-1, false, false>
    BaseDerivType;

  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  PowFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      df ( BaseDerivType((FloatType)p, PowM1Type(x)) ),
      value (df.right.value * x.value) {}
  template<bool U = refArg>
  PowFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
      df ( BaseDerivType((FloatType)p, PowM1Type(x)) ),
      value (df.right.value * x.value) {}

  const BaseDerivType df;
  const FloatType value;  

  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) * df )
  {
    return getD<derivOrder>(x, index) * df;
  }
  
  typedef decltype(((PowFad*)NULL)->d(0)) DerivType;
};

// -Fad with derivOrder > 1, p == 1
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg>
struct PowFad <FloatType, ArgType, derivOrder, 1, refArg>
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;
    
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  PowFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      value (x.value) {}
  template<bool U = refArg>
  PowFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
     value (x.value) {}
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(x, index))
  {
    return getD<derivOrder>(x, index);
  }

  typedef decltype(((PowFad*)NULL)->d(0)) DerivType;
  
  const FloatType & value;
};

// -Fad  with derivOrder == 1, p > 1
template <typename FloatType, typename ArgType, size_t p, bool refArg>
struct PowFad<FloatType, ArgType, 1, p, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef a short name for the derivatives' type
  typedef FloatType DerivType;

  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  PowFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      xPM1 (power<p-1>(x.value)),
      df ((FloatType)p * xPM1),
      value (x.value * xPM1) {}
  template<bool U = refArg>
  PowFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
     xPM1 (power<p-1>(x.value)),
     df ((FloatType)p * xPM1),
     value (x.value * xPM1) {}
  
  // return the derivative as a FloatType
  FloatType d(size_t index) const {
    return x.hasD(index) ? getD<1>(x, index) * df : (FloatType)0.0;
  }
  
  // return the derivative as a FloatType
  FloatType D(size_t index) const { return getD<1>(x, index) * df; }
    
  const FloatType xPM1;
  const FloatType df;
  const FloatType value;
};

// -Fad  with derivOrder == 1, p == 1
template <typename FloatType, typename ArgType, bool refArg>
struct PowFad<FloatType, ArgType, 1, 1, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef a short name for the derivatives' type
  typedef FloatType DerivType;
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  PowFad(ArgType && op, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(op)),
      value (x.value) {}
  template<bool U = refArg>
  PowFad(const ArgType & op, typename std::enable_if<U>::type* = 0)
   : OpType (op),
     value (x.value) {}
  
  // return the derivative as a FloatType
  FloatType d(size_t index) const {
    return x.hasD(index) ? getD<1>(x, index) : (FloatType)0.0;
  }
  
  // return the derivative as a FloatType
  FloatType D(size_t index) const { return getD<1>(x, index); }
  
  const FloatType & value;
};



template <size_t p, FAD_TYPENAME(T)>
inline
PowFad<typename T::Float_T, T, T::getDerivOrder(), p, false>
power(T && x)
{
  return PowFad<typename T::Float_T, T, T::getDerivOrder(), p, false>
    (std::move(x));
}
template <size_t p, FAD_TYPENAME(T)>
inline
PowFad<typename T::Float_T, T, T::getDerivOrder(), p, true>
power(const T & x)
{
  return PowFad<typename T::Float_T, T, T::getDerivOrder(), p, true> (x);
}



//////////////////////////// DivFad - Division class //////////////////////////

// Fad / Fad with derivOrder > 1
template <typename FloatType, typename LeftType, typename RightType,
           Order_T derivOrder, bool refLeft=false, bool refRight=false>
struct DivFad
  : public BinaryOpFad<FloatType, LeftType, RightType, derivOrder,
                       refLeft, refRight>
{
  // compute derivative according to formula
  // d/dt (left/right) = (d/dt left - (left/right) * d/dt right) / right
  
  // typedef the inner division
  typedef DivFad<FloatType, LeftType, RightType, derivOrder-1, true, true>
    InnerDivType;
  
  // typedef the base type and allow local access to inherited member variables
  typedef BinaryOpFad<FloatType, LeftType, RightType, derivOrder,
                      refLeft, refRight> OpType;
  using OpType::left;  using OpType::right;
  
  // define the constructor
  template<bool L = refLeft, bool R = refRight>
  DivFad(LeftType && op1, RightType && op2,
         typename std::enable_if<!L&!R>::type* = 0)
    : OpType(std::move(op1), std::move(op2)),
      div(left, right),
      value (div.value) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(LeftType && op1, const RightType & op2,
         typename std::enable_if<!L&R>::type* = 0)
    : OpType(std::move(op1), op2),
      div(left, right),
      value (div.value) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(const LeftType & op1, RightType && op2,
        typename std::enable_if<L&!R>::type* = 0)
    : OpType(op1, std::move(op2)),
      div(left, right),
      value (div.value) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(const LeftType & op1, const RightType & op2,
         typename std::enable_if<L&R>::type* = 0)
    : OpType(op1, op2),
      div(left, right),
      value (div.value) {}

  InnerDivType div;
  const FloatType & value;
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype((getD<derivOrder>(left, index)
                 - div * getD<derivOrder>(right, index)) / right)
  {
    return (getD<derivOrder>(left, index)
            - div * getD<derivOrder>(right, index)) / right;
  }
  
  typedef decltype(((DivFad*)NULL)->d(0)) DerivType;
  
};

// Fad / Float with derivOrder > 1
template <typename FloatType, typename LeftType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct DivFad<FloatType, LeftType, FloatType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, derivOrder,
                       refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(DivFad, LeftType, FloatType, derivOrder,
                     left.value / right)
  
  // return the derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype(getD<derivOrder>(left, index) / right)
  {
    return getD<derivOrder>(left, index) / right;
  }
  
  typedef decltype(((DivFad*)NULL)->d(0)) DerivType;
};

// Float / Fad with derivOrder > 1
template <typename FloatType, typename RightType, Order_T derivOrder,
          bool refLeft, bool refRight>
struct DivFad<FloatType, FloatType, RightType, derivOrder, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, derivOrder,
                       refLeft, refRight>
{
  // compute derivative according to formula
  // d/dt (left/right) = ((-left)/power<2>(right)) * d/dt right
  // d/dt (left/right) = -(left/right) / right * d/dt(right)
  
  // typedef the base type and allow local access to inherited member variables
  typedef BinaryOpFad<FloatType, FloatType, RightType, derivOrder,
                      refLeft, refRight> OpType;
  using OpType::left;  using OpType::right;

  const FloatType value;

  // typedef the derivative formula
  typedef decltype((-left)/power<2>(wrap<derivOrder-1>(right))) BaseDerivType;
  BaseDerivType df;

  // define the constructor
  template<bool L = refLeft, bool R = refRight>
  DivFad(FloatType && op1, RightType && op2,
         typename std::enable_if<!L&!R>::type* = 0)
    : OpType(std::move(op1), std::move(op2)),
      value ( left / right.value ),
      df ( (-left) / power<2>(wrap<derivOrder-1>(right)) ) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(FloatType && op1, const RightType & op2,
         typename std::enable_if<!L&R>::type* = 0)
    : OpType(std::move(op1), op2),
      value ( left / right.value ),
      df ( (-left) / power<2>(wrap<derivOrder-1>(right)) ) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(const FloatType & op1, RightType && op2,
        typename std::enable_if<L&!R>::type* = 0)
    : OpType(op1, std::move(op2)),
      value ( left / right.value ),
      df ( (-left) / power<2>(wrap<derivOrder-1>(right)) ) {}
  template<bool L = refLeft, bool R = refRight>
  DivFad(const FloatType & op1, const RightType & op2,
         typename std::enable_if<L&R>::type* = 0)
    : OpType(op1, op2),
      value ( left / right.value ),
      df ( (-left) / power<2>(wrap<derivOrder-1>(right)) ) {}
      
  // return the derivative as a Fad-type
  auto d(size_t index) const -> decltype(getD<derivOrder>(right, index) * df)
  {
    return getD<derivOrder>(right, index) * df;
  }
  
  typedef decltype(((DivFad*)NULL)->d(0)) DerivType;
};

// Fad / Fad with derivOrder == 1
template <typename FloatType, typename LeftType, typename RightType,
          bool refLeft, bool refRight>
struct DivFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(DivFad, LeftType, RightType, 1, left.value / right.value)

  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index)
      ? right.hasD(index) ? (getD<1>(left, index)
                             - value * getD<1>(right, index)) / right.value
                          : getD<1>(left, index) / right.value
      : right.hasD(index) ? -value * getD<1>(right, index) / right.value
                          : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return left.hasD(index)
      ? right.hasD(index) ? (getD<1>(left, index)
                             - value * getD<1>(right, index)) / right.value
                          : getD<1>(left, index) / right.value
      : -value * getD<1>(right, index) / right.value;
  }
  
  typedef FloatType DerivType;
};

// Fad / Float with derivOrder == 1
template <typename FloatType, typename LeftType, bool refLeft, bool refRight>
struct DivFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, LeftType, FloatType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(DivFad, LeftType, FloatType, 1, left.value / right)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return left.hasD(index) ? getD<1>(left, index) / right : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return getD<1>(left, index) / right;
  }
  
  typedef FloatType DerivType;
};

// Float / Fad with derivOrder == 1
template <typename FloatType, typename RightType, bool refLeft, bool refRight>
struct DivFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
  : public BinaryOpFad<FloatType, FloatType, RightType, 1, refLeft, refRight>
{
  // typedef the base type,
  // allow local access to inherited member variables,
  // and define constructor
  BINARY_CONSTRUCTOR(DivFad, FloatType, RightType, 1, left / right.value)
  
  // return the derivative as a FloatType
  FloatType d(const size_t & index) const
  {
    return right.hasD(index) ? -value * getD<1>(right, index) / right.value
                             : (FloatType)0.0;
  }
  FloatType D(const size_t & index) const
  {
    return -value * getD<1>(right, index) / right.value;
  }
  
  typedef FloatType DerivType;
};


DEFINE_FAD_ARITHMETIC_OPERATOR(Div, / )

#endif
