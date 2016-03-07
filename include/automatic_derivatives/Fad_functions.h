#ifndef Fad_h
#error Do not include this file directly, instead include Fad.h
#endif
#ifndef Fad_functions_h
#define Fad_functions_h


#include <cmath>



//////////// Teach Fad how to differentiate basic special functions////////////



// define a macro that's useful for declaring helper Fad classes and defining
//   functions that act on Fad-types
#define DECLARE_FUNC(funcName) \
template <typename FloatType, typename ArgType, Order_T derivOrder, bool refArg> \
struct funcName##Fad; \
template <typename XType, \
          typename = typename std::enable_if<XType::isFad>::type> \
inline funcName##Fad<typename XType::Float_T, XType, XType::getDerivOrder(), \
                     false> \
funcName(XType && x) \
{ \
  return funcName##Fad<typename XType::Float_T, XType, XType::getDerivOrder(), \
                       false> (std::move(x)); \
} \
template <typename XType, \
          typename = typename std::enable_if<XType::isFad>::type> \
inline funcName##Fad<typename XType::Float_T, XType, XType::getDerivOrder(), \
                     true> \
funcName(const XType & x) \
{ \
  return funcName##Fad<typename XType::Float_T, XType, XType::getDerivOrder(), \
                       true> (x); \
}

// list of functions that are supported in this file:

// These functions have tricks that enable them to be computed more quickly
// (the derivative is related to the value of the original function)
DECLARE_FUNC(exp)
DECLARE_FUNC(log)
DECLARE_FUNC(sqrt)

// sin, cos, sinh, cosh use tricks with special classes TrigFad, HTrigFad
// They are not declared here because they don't fit the usual pattern
// DECLARE_FUNC(sin)
// DECLARE_FUNC(cos)
// DECLARE_FUNC(sinh)
// DECLARE_FUNC(cosh)

// These functions all have derivatives of form df(x)/dt = dx/dt / g(x)
// They are computed with a special class
DECLARE_FUNC(log10)
DECLARE_FUNC(tan)
DECLARE_FUNC(tanh)
DECLARE_FUNC(atan)
DECLARE_FUNC(atanh)
DECLARE_FUNC(asin)
DECLARE_FUNC(acos)
DECLARE_FUNC(asinh)
DECLARE_FUNC(acosh)

// pow(x,p) is defined as exp(log(x) * p)



//////////////////// define expFad for computing exp(x) ///////////////////////
// d/dx exp(x) = exp(x)
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct expFad : UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{    
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;
 
  // typedef the base derivative type and declare df
  typedef expFad<FloatType, ArgType, derivOrder-1, true>
    BaseDerivType;
  const BaseDerivType df;
  
  // declare the value
  const FloatType & value;
 
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  expFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(x) ),
      value (df.value) {}
  template<bool U = refArg>
  expFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
   : OpType (arg),
      df ( BaseDerivType(x) ),
      value (df.value) {}

  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) * df )
  {
    return getD<derivOrder>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef decltype(((expFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, bool refArg>
struct expFad<FloatType, ArgType, 1, refArg>
  : UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // declare the value
  const FloatType value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  expFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value ( exp(x.value) ) {}
  template<bool U = refArg>
  expFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
   : OpType (arg),
     value ( exp(x.value) ) {}
  
  // return derivative as FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) * value : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) * value;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};



//////////////////// define logFad for computing log(x) ///////////////////////
// d/dx log(x) = 1/x
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct logFad : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;

  // declare the value
  const FloatType value;

  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  logFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value (log(x.value)) {}
  template<bool U = refArg>
  logFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      value (log(x.value)) {}

  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) / x )
  {
    return getD<derivOrder>(x, index) / x;
  }
  
  // typedef the DerivType
  typedef decltype(((logFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, bool refArg>
struct logFad<FloatType, ArgType, 1, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // declare the value
  const FloatType value;

  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  logFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value (log(x.value)) {}
  template<bool U = refArg>
  logFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      value (log(x.value)) {}

  // return derivative as FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) / x.value : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) / x.value;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};



/////////////////// define sqrtFad for computing sqrt(x) //////////////////////
// d/dx sqrt(x) = 1/(sqrt(x) * 2.0)
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct sqrtFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
 // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;
  
  // typedef the base derivative type and declare df
  typedef sqrtFad<FloatType, ArgType, derivOrder-1, true> RootType;
  typedef MulFad<FloatType, RootType, FloatType, derivOrder-1,
                 false, true> BaseDerivType;
  const BaseDerivType df;
  
  // declare value
  const FloatType & value;
 
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  sqrtFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(RootType(x), (FloatType)2.0) ),
      value ( df.left.value ) {}
  template<bool U = refArg>
  sqrtFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
   : OpType (arg),
      df ( BaseDerivType(RootType(x), (FloatType)2.0) ),
      value ( df.left.value ) {}

  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) / df )
  {
    return getD<derivOrder>(x, index) / df;
  }
  
  // typedef the DerivType
  typedef decltype(((sqrtFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, bool refArg>
struct sqrtFad<FloatType, ArgType, 1, refArg>
 : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // declare value
  const FloatType value;
  
  // declare base derivative
  const FloatType df;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  sqrtFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value ( sqrt(x.value) ),
      df ( 2.0 * value ) {}
  template<bool U = refArg>
  sqrtFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
   : OpType (arg),
      value ( sqrt(x.value) ),
      df ( 2.0 * value ) {}
  
  // return derivative as FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) / df : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) / df;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};



/////////// Helper functions/classes for computing sin(x) and cos(x) //////////

enum class TrigPhase {sin, cos, msin, mcos};
template <TrigPhase phase>
constexpr TrigPhase
nextTrigPhase(void)
{
  return TrigPhase::cos;
}
template <>
constexpr TrigPhase
nextTrigPhase<TrigPhase::cos>(void)
{
  return TrigPhase::msin;
}
template <>
constexpr TrigPhase
nextTrigPhase<TrigPhase::msin>(void)
{
  return TrigPhase::mcos;
}
template <>
constexpr TrigPhase
nextTrigPhase<TrigPhase::mcos>(void)
{
  return TrigPhase::sin;
}



template <typename FloatType, TrigPhase phase>
class trigFunc
{
  public:
    inline FloatType f(FloatType x) { return sin(x); }
};
template <typename FloatType>
class trigFunc<FloatType, TrigPhase::cos>
{
  public:
    inline FloatType f(FloatType x) { return cos(x); }
};
template <typename FloatType>
class trigFunc<FloatType, TrigPhase::msin>
{
  public:
    inline FloatType f(FloatType x) { return -sin(x); }
};
template <typename FloatType>
class trigFunc<FloatType, TrigPhase::mcos>
{
  public:
    inline FloatType f(FloatType x) { return -cos(x); }
};



/////////////// define TrigFad for computing sin(x) and cos(x) ////////////////
// d/dx sin(x) = cos(x), d/dx cos(x) = -sin(x)
template <typename FloatType, typename ArgType, Order_T derivOrder,
          TrigPhase phase, bool refArg=false>
struct TrigFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;

  // typedef the base derivative type and declare df
  typedef TrigFad<FloatType, ArgType, derivOrder-1, nextTrigPhase<phase>(),
                  true> BaseDerivType;
  const BaseDerivType df;
  
  // declare the value
  const FloatType value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  TrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(x) ),
      value ( -df.df.value ) {}
  template<bool U = refArg>
  TrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      df ( BaseDerivType(x) ),
      value ( -df.df.value ) {}
  
  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) * df )
  {
    return getD<derivOrder>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef decltype(((TrigFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, TrigPhase phase, bool refArg>
struct TrigFad<FloatType, ArgType, 2, phase, refArg>
  : public UnaryOpFad<FloatType, ArgType, 2, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 2, refArg> OpType;
  using OpType::x;

  // typedef the base derivative type and declare df
  typedef TrigFad<FloatType, ArgType, 1, nextTrigPhase<phase>(), true>
    BaseDerivType;
  const BaseDerivType df;
  
  // declare the value
  const FloatType value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  TrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(x) ),
      value ( -df.df ) {}
  template<bool U = refArg>
  TrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      df ( BaseDerivType(x) ),
      value ( -df.df ) {}
  
  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<2>(x, index) * df )
  {
    return getD<2>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef decltype(((TrigFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, TrigPhase phase, bool refArg>
struct TrigFad<FloatType, ArgType, 1, phase, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // declare the value
  const FloatType value;
  
  // declare the base derivative
  const FloatType df;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  TrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value ( trigFunc<FloatType, phase>().f(x.value) ),
      df ( trigFunc<FloatType, nextTrigPhase<phase>()>().f(x.value) ) {}
  template<bool U = refArg>
  TrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      value ( trigFunc<FloatType, phase>().f(x.value) ),
      df ( trigFunc<FloatType, nextTrigPhase<phase>()>().f(x.value) ) {}
  
  // return derivative as FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) * df : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};



template <FAD_TYPENAME(XType)>
inline TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
               TrigPhase::sin, false>
sin(XType && x)
{
  return TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                 TrigPhase::sin, false> (std::move(x));
}
template <FAD_TYPENAME(XType)>
inline TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
               TrigPhase::sin, true>
sin(const XType & x)
{
  return TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                 TrigPhase::sin, true> (x);
}
template <FAD_TYPENAME(XType)>
inline TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
               TrigPhase::cos, false>
cos(XType && x)
{
  return TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                 TrigPhase::cos, false> (std::move(x));
}
template <FAD_TYPENAME(XType)>
inline TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
               TrigPhase::cos, true>
cos(const XType & x)
{
  return TrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                 TrigPhase::cos, true> (x);
}


////////// Helper functions/classes for computing sinh(x) and cosh(x) /////////

enum class HTrigPhase {sinh, cosh};
template <HTrigPhase phase>
constexpr HTrigPhase
nextHTrigPhase(void)
{
  return HTrigPhase::cosh;
}
template <>
constexpr HTrigPhase
nextHTrigPhase<HTrigPhase::cosh>(void)
{
  return HTrigPhase::sinh;
}



template <typename FloatType, HTrigPhase phase>
class hTrigFunc
{
  public:
    inline FloatType f(FloatType x) { return sinh(x); }
};
template <typename FloatType>
class hTrigFunc<FloatType, HTrigPhase::cosh>
{
  public:
    inline FloatType f(FloatType x) { return cosh(x); }
};



////////////// define HTrigFad for computing sinh(x) and cosh(x) //////////////
// d/dx sinh(x) = cosh(x), d/dx cosh(x) = sinh(x)
template <typename FloatType, typename ArgType, Order_T derivOrder,
          HTrigPhase phase, bool refArg=false>
struct HTrigFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;

  // typedef the base derivative type and declare df
  typedef HTrigFad<FloatType, ArgType, derivOrder-1, nextHTrigPhase<phase>(),
                   true> BaseDerivType;
  const BaseDerivType df;
  
  // declare the value
  const FloatType & value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  HTrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(x) ),
      value ( df.df.value ) {}
  template<bool U = refArg>
  HTrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      df ( BaseDerivType(x) ),
      value ( df.df.value ) {}
  
  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) * df )
  {
    return getD<derivOrder>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef decltype(((HTrigFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, HTrigPhase phase, bool refArg>
struct HTrigFad<FloatType, ArgType, 2, phase, refArg>
  : public UnaryOpFad<FloatType, ArgType, 2, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 2, refArg> OpType;
  using OpType::x;

  // typedef the base derivative type and declare df
  typedef HTrigFad<FloatType, ArgType, 1, nextHTrigPhase<phase>(),
                   true> BaseDerivType;
  const BaseDerivType df;
  
  // declare the value
  const FloatType & value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  HTrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( BaseDerivType(x) ),
      value ( df.df ) {}
  template<bool U = refArg>
  HTrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      df ( BaseDerivType(x) ),
      value ( df.df ) {}
  
  // return derivative as a Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<2>(x, index) * df )
  {
    return getD<2>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef decltype(((HTrigFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, HTrigPhase phase, bool refArg>
struct HTrigFad<FloatType, ArgType, 1, phase, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;

  // declare the value
  const FloatType value;
  
  // declare the base derivative
  const FloatType df;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  HTrigFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      value ( hTrigFunc<FloatType, phase>().f(x.value) ),
      df ( hTrigFunc<FloatType, nextHTrigPhase<phase>()>().f(x.value) ) {}
  template<bool U = refArg>
  HTrigFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      value ( hTrigFunc<FloatType, phase>().f(x.value) ),
      df ( hTrigFunc<FloatType, nextHTrigPhase<phase>()>().f(x.value) ) {}
  
  // return derivative as FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) * df : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) * df;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};


template <FAD_TYPENAME(XType)>
inline HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                HTrigPhase::sinh, false>
sinh(XType && x)
{
  return HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                  HTrigPhase::sinh, false> (std::move(x));
}
template <FAD_TYPENAME(XType)>
inline HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                HTrigPhase::sinh, true>
sinh(const XType & x)
{
  return HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                  HTrigPhase::sinh, true> (x);
}

template <FAD_TYPENAME(XType)>
inline HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                HTrigPhase::cosh, false>
cosh(XType && x)
{
  return HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                  HTrigPhase::cosh, false> (std::move(x));
}
template <FAD_TYPENAME(XType)>
inline HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                HTrigPhase::cosh, true>
cosh(const XType & x)
{
  return HTrigFad<typename XType::Float_T, XType, XType::getDerivOrder(),
                  HTrigPhase::cosh, true> (x);
}



//// define struct for computing derivatives of form df(x)/dt = dx/dt / g(x) ///

template <typename FloatType, typename ArgType, Order_T derivOrder,
          typename EvalFunctor, bool refValue=false, bool refArg=false>
struct FuncUnaryOpFad
  : public UnaryOpFad<FloatType, ArgType, derivOrder, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, derivOrder, refArg> OpType;
  using OpType::x;

  // save a version of x with lower derivOrder
  typedef decltype(wrap<derivOrder-1>(x)) WrapXType;
  const WrapXType wrapX;
  
  // typedef the BaseDerivType
  typedef decltype(EvalFunctor::df(wrapX)) BaseDerivType;
  // declare the base derivative
  const BaseDerivType df;

  // declare the value
  typename RefArgType<FloatType, refValue>::type value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  FuncUnaryOpFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      wrapX ( wrap<derivOrder-1>(x) ),
      df ( EvalFunctor::df(wrapX) ),
      value ( EvalFunctor::f(x.value) ) {}
  template<bool U = refArg>
  FuncUnaryOpFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      wrapX ( wrap<derivOrder-1>(x) ),
      df ( EvalFunctor::df(wrapX) ),
      value ( EvalFunctor::f(x.value) ) {}

  // return derivative as Fad-type
  auto d(const size_t & index) const
    -> decltype( getD<derivOrder>(x, index) / df )
  {
    return getD<derivOrder>(x, index) / df;
  }
  
  // typedef the DerivType
  typedef decltype(((FuncUnaryOpFad*)NULL)->d(0)) DerivType;
};

template <typename FloatType, typename ArgType, typename EvalFunctor,\
          bool refValue, bool refArg>
struct FuncUnaryOpFad<FloatType, ArgType, 1, EvalFunctor, refValue, refArg>
  : public UnaryOpFad<FloatType, ArgType, 1, refArg>
{
  // typedef the base type and allow local access to inherited member variables
  typedef UnaryOpFad<FloatType, ArgType, 1, refArg> OpType;
  using OpType::x;
  
  // declare the base derivative
  const FloatType df;

  // declare the value
  typename RefArgType<FloatType, refValue>::type value;
  
  // define constructor, with and without pass-by-reference
  template<bool U = refArg>
  FuncUnaryOpFad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)),
      df ( EvalFunctor::df(x.value) ),
      value ( EvalFunctor::f(x.value) ) {}
  template<bool U = refArg>
  FuncUnaryOpFad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg),
      df ( EvalFunctor::df(x.value) ),
      value ( EvalFunctor::f(x.value) ) {}
  
  // return derivative as a FloatType
  FloatType d(const size_t & index) const {
    return x.hasIndex(index) ? getD<1>(x, index) / df : (FloatType)0.0;
  }

  FloatType D(const size_t & index) const {
    return getD<1>(x, index) / df;
  }
  
  // typedef the DerivType
  typedef FloatType DerivType;
};


//////// Define macro to pass derivative information to FuncUnaryOpFad ////////

#define DEFINE_DERIV(funcName, dFunc, invDeriv) \
struct funcName##EvalFunctor \
{ \
  template <typename T> \
  static auto f(const T & x) -> decltype( funcName (x) ) \
  { return funcName (x); } \
  template <typename T> \
  static auto df(const T & x) -> decltype( dFunc ) { return dFunc; } \
}; \
template <typename FloatType, typename ArgType, Order_T derivOrder, \
          bool refArg=false> \
struct funcName##Fad \
  : public FuncUnaryOpFad<FloatType, ArgType, derivOrder, \
                          funcName##EvalFunctor, false, refArg> \
{ \
  typedef FuncUnaryOpFad<FloatType, ArgType, derivOrder, \
                         funcName##EvalFunctor, false, refArg> \
    OpType; \
  template<bool U = refArg> \
  funcName##Fad(ArgType && arg, typename std::enable_if<!U>::type* = 0) \
    : OpType (std::move(arg)) {} \
  template<bool U = refArg> \
  funcName##Fad(const ArgType & arg, typename std::enable_if<U>::type* = 0) \
    : OpType (arg) {} \
};


/////////////////////////// DERIVATIVE DEFINITIONS ////////////////////////////
/*
struct log10EvalFunctor
{
  template <typename T>
  static auto f(const T & x) -> decltype( log10(x) )
  { return log10(x); }

  template <typename T>
  static auto df(const T & x) -> decltype( x * log(10.0) )
  { return x * log(10.0); }
};
template <typename FloatType, typename ArgType, Order_T derivOrder,
          bool refArg=false>
struct log10Fad \
  : public FuncUnaryOpFad<FloatType, ArgType, derivOrder,
                          log10EvalFunctor, true, false, refArg> \
{
  typedef FuncUnaryOpFad<FloatType, ArgType, derivOrder,
                         log10EvalFunctor, true, false, refArg>
    OpType;
  template<bool U = refArg>
  log10Fad(ArgType && arg, typename std::enable_if<!U>::type* = 0)
    : OpType (std::move(arg)) {}
  template<bool U = refArg>
  log10Fad(const ArgType & arg, typename std::enable_if<U>::type* = 0)
    : OpType (arg) {}
};
*/


// d/dx log10(x) = 1.0 / (x * log(10.0))
DEFINE_DERIV(log10, x * log(10.0), true);
// d/dx tan(x) = 1.0 / power<2>(cos(x))
DEFINE_DERIV(tan, power<2>(cos(x)), true);
// d/dx tanh(x) = 1.0 / power<2>(cosh(x))
DEFINE_DERIV(tanh, power<2>(cosh(x)), true);
// d/dx atan(x) = 1 / (1 + power<2>(x))
DEFINE_DERIV(atan, 1.0 + power<2>(x), true);
// d/dx atanh(x) = 1 / (1 - power<2>(x))
DEFINE_DERIV(atanh, 1.0 - power<2>(x), true);
// d/dx asin(x) = 1.0 / sqrt(1.0 - power<2>(x))
DEFINE_DERIV(asin, sqrt(1.0 - power<2>(x)), true);
// d/dx acos(x) = -1.0 / sqrt(1.0 - power<2>(x))
DEFINE_DERIV(acos, -sqrt(1.0 - power<2>(x)), true);
// d/dx asinh(x) = 1.0 / sqrt(power<2>(x) + 1.0)
DEFINE_DERIV(asinh, sqrt(1.0 + power<2>(x)), true);
// d/dx acosh(x) = 1.0 / sqrt(power<2>(x) - 1.0)
DEFINE_DERIV(acosh, sqrt(power<2>(x) - 1.0), true);



/////////////////////////////// pow functions /////////////////////////////////

template <FAD_TYPENAME(BaseType), FAD_TYPENAME(ExpType)>
inline auto
pow(BaseType && x, ExpType && p)
 -> decltype(exp(log(std::move(x)) * std::move(p)))
{
  static_assert(std::is_same<typename BaseType::Float_T,
                             typename ExpType::Float_T>::value,
                "x and p have incompatible float types");
  static_assert(BaseType::getDerivOrder() == ExpType::getDerivOrder(),
                "x and p have incompatible derivative orders");

  return exp(log(std::move(x)) * std::move(p));
}
template <FAD_TYPENAME(BaseType), FAD_TYPENAME(ExpType)>
inline auto
pow(BaseType && x, const ExpType & p)
 -> decltype(exp(log(std::move(x)) * p))
{
  static_assert(std::is_same<typename BaseType::Float_T,
                             typename ExpType::Float_T>::value,
                "x and p have incompatible float types");
  static_assert(BaseType::getDerivOrder() == ExpType::getDerivOrder(),
                "x and p have incompatible derivative orders");

  return exp(log(std::move(x)) * p);
}
template <FAD_TYPENAME(BaseType), FAD_TYPENAME(ExpType)>
inline auto
pow(const BaseType & x, ExpType && p)
 -> decltype(exp(log(x) * std::move(p)))
{
  static_assert(std::is_same<typename BaseType::Float_T,
                             typename ExpType::Float_T>::value,
                "x and p have incompatible float types");
  static_assert(BaseType::getDerivOrder() == ExpType::getDerivOrder(),
                "x and p have incompatible derivative orders");

  return exp(log(x) * std::move(p));
}
template <FAD_TYPENAME(BaseType), FAD_TYPENAME(ExpType)>
inline auto
pow(const BaseType & x, const ExpType & p)
 -> decltype(exp(log(x) * p))
{
  static_assert(std::is_same<typename BaseType::Float_T,
                             typename ExpType::Float_T>::value,
                "x and p have incompatible float types");
  static_assert(BaseType::getDerivOrder() == ExpType::getDerivOrder(),
                "x and p have incompatible derivative orders");

  return exp(log(x) * p);
}


template <FAD_TYPENAME(BaseType), FLOAT_TYPENAME(FloatType)>
inline auto
pow(BaseType && x, FloatType p)
 -> decltype(exp(log(std::move(x)) * (typename BaseType::Float_T)p))
{
  return exp(log(std::move(x)) * (typename BaseType::Float_T)p);
}

template <FAD_TYPENAME(BaseType), FLOAT_TYPENAME(FloatType)>
inline auto
pow(const BaseType & x, FloatType p)
 -> decltype(exp(log(x) * (typename BaseType::Float_T)p))
{
  return exp(log(x) * (typename BaseType::Float_T)p);
}


template <FLOAT_TYPENAME(FloatType), FAD_TYPENAME(ExpType)>
inline auto
pow(FloatType x, ExpType && p)
 -> decltype(exp(log((typename ExpType::Float_T)x) * std::move(p)))
{
  return exp(log((typename ExpType::Float_T)x) * std::move(p));
}

template <FLOAT_TYPENAME(FloatType), FAD_TYPENAME(ExpType)>
inline auto
pow(FloatType x, const ExpType & p)
 -> decltype(exp(log((typename ExpType::Float_T)x) * p))
{
  return exp(log((typename ExpType::Float_T)x) * p);
}

#endif
