// -*- C++ -*-
// $Id: diagnostic.t,v 1.2 2003/12/11 10:40:10 hs Exp $


template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::open(const std::string &fname)
{
  output.open(fname.c_str());
  output.precision(14);
}

template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::write()
{
   output << *field;
}

template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::close()
{
  output.close();
}

template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::setField(Type *fld)
{
  field = fld;
}


