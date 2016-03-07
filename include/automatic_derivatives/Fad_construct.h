#ifndef Fad_h
#error Do not include this file directly, instead include Fad.h
#endif
#ifndef Fad_construct_h
#define Fad_construct_h


/////////////////////////////// Helper Functions //////////////////////////////


////////////////////////////// Fad Constructors ///////////////////////////////

template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder>::Fad()
{
  // default constructor
}


template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder>::Fad(const FloatType & assignValue)
  : value (assignValue)
{
  // construct and assign value
}


template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder>::Fad(const FloatType & assignValue,
                                const size_t & derivIndex)
  : value (assignValue), 
    derivatives ( {{derivIndex, DerivType(1.0)}} )
{
  // construct, assign value and set derivative index
}



template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder>::Fad(const Fad<FloatType, derivOrder> & other)
  : value (other.value), 
    derivatives (other.derivatives)
{
  // copy constructor
}



template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder>::Fad(Fad<FloatType, derivOrder> && other)
  : value (std::move(other.value))
{
  // r-value constructor
  std::swap(derivatives, other.derivatives);
}



template <typename FloatType>
Fad<FloatType, 1>::Fad()
{
  // default constructor
}



template <typename FloatType>
Fad<FloatType, 1>::Fad(const FloatType & assignValue)
  : value (assignValue)
{
  // construct and assign value
}



template <typename FloatType>
Fad<FloatType, 1>::Fad(const FloatType & assignValue,
                       const size_t & derivIndex)
  : value (assignValue), 
    derivatives ( {{derivIndex, (FloatType)1.0}} )
{
  // construct, assign value and set derivative index
}



template <typename FloatType>
Fad<FloatType, 1>::Fad(const Fad<FloatType, 1> & other)
  : value (other.value), 
    derivatives (other.derivatives)
{
  // copy constructor
}



template <typename FloatType>
Fad<FloatType, 1>::Fad(Fad<FloatType, 1> && other)
  : value (std::move(other.value))
{
  // r-value copy constructor

  std::swap(derivatives, other.derivatives);
}



//////////////////////// Fad Allocation and Assignment ////////////////////////



template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder> &
Fad<FloatType, derivOrder>::operator=(const Fad<FloatType, derivOrder> & other)
{
  // assign other FadType to *this
  
  // copy value
  value = other.value;
  
  derivatives = other.derivatives;
  
  return *this;
}



template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder> &
Fad<FloatType, derivOrder>::operator=(Fad<FloatType, derivOrder> && other)
{
  // assign temporary other FadType to *this
  
  // move value
  value = std::move(other.value);
  
  // swap derivatives, since temporary will be destroyed
  std::swap(derivatives, other.derivatives);
  
  return *this;
}



template <typename FloatType, unsigned char derivOrder>
Fad<FloatType, derivOrder> &
Fad<FloatType, derivOrder>::operator=(const FloatType & other)
{
  // assign other FloatType to this->value
  value = other;
  
  return *this;
}



template <typename FloatType>
Fad<FloatType, 1> &
Fad<FloatType, 1>::operator=(const Fad<FloatType, 1> & other)
{
  // assign other FadType to *this
  
  // copy value
  value = other.value;
  
  derivatives = other.derivatives;
  
  return *this;
}



template <typename FloatType>
Fad<FloatType, 1> &
Fad<FloatType, 1>::operator=(Fad<FloatType, 1> && other)
{
  // assign temporary other FadType to *this
  
  // move value
  value = std::move(other.value);
  
  // swap derivatives, since temporary will be destroyed
  std::swap(derivatives, other.derivatives);
  
  return *this;
}



template <typename FloatType>
Fad<FloatType, 1> &
Fad<FloatType, 1>::operator=(const FloatType & other)
{
  // assign other FloatType to this->value
  value = other;
  
  return *this;
}



//////////////////////// Fad Miscellaneous Functions //////////////////////////

template <typename FloatType, unsigned char derivOrder>
void
Fad<FloatType, derivOrder>::clear(void)
{
  // clear the derivative information
  derivatives.clear();
}
  

template <typename FloatType, unsigned char derivOrder>
void
Fad<FloatType, derivOrder>::setDerivativeIndex(const size_t & derivIndex)
{
  // set the derivative index for this variable
  derivatives.clear();
  derivatives[derivIndex] = DerivType(1.0);
}


template <typename FloatType, unsigned char derivOrder>
bool
Fad<FloatType, derivOrder>::hasD(const size_t & index) const
{
  // check if an index has a nonzero derivative
  return derivatives.count(index);
}


template <typename FloatType, unsigned char derivOrder>
const typename Fad<FloatType, derivOrder>::DerivType &
Fad<FloatType, derivOrder>::D(const size_t & index) const
{
  // return the derivative of the Fad without checking indices first
  return derivatives.at(index);
}


template <typename FloatType, unsigned char derivOrder>
typename Fad<FloatType, derivOrder>::DerivType
Fad<FloatType, derivOrder>::d(const size_t & index) const {
  // return the derivative of the Fad as another Fad object
  if(derivatives.count(index)) {
    return derivatives.at(index);
  }
  else {
    return DerivType(0.0);
  }
}


template <typename FloatType, unsigned char derivOrder>
std::set<size_t>
Fad<FloatType, derivOrder>::getIndices(void) const {
  // return the indices to non-zero derivatives
  std::set<size_t> indices;
  for(const auto mapObj : derivatives) {
    indices.insert(mapObj.first);
  }
  return indices;
}



template <typename FloatType>
void
Fad<FloatType, 1>::clear(void)
{
  // clear the derivative information
  derivatives.clear();
}
  

template <typename FloatType>
void
Fad<FloatType, 1>::setDerivativeIndex(const size_t & derivIndex)
{
  // set the derivative index for this variable
  derivatives.clear();
  derivatives[derivIndex] = DerivType(1.0);
}


template <typename FloatType>
bool
Fad<FloatType, 1>::hasD(const size_t & index) const
{
  // check if an index has a nonzero derivative
  return derivatives.count(index);
}


template <typename FloatType>
const FloatType &
Fad<FloatType, 1>::D(const size_t & index) const
{
  // return the derivative of the Fad without checking indices first
  return derivatives.at(index);
}


template <typename FloatType>
FloatType
Fad<FloatType, 1>::d(const size_t & index) const {
  // return the derivative of the Fad as another Fad object
  return derivatives.count(index) ?  derivatives.at(index) : (FloatType)0.0;
}


template <typename FloatType>
std::set<size_t>
Fad<FloatType, 1>::getIndices(void) const {
  // return the indices to non-zero derivatives
  std::set<size_t> indices;
  for(const auto mapObj : derivatives) {
    indices.insert(mapObj.first);
  }
  return indices;
}
#endif
