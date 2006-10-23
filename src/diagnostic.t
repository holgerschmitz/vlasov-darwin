// -*- C++ -*-
// $Id: diagnostic.t,v 1.1.1.1 2006/05/11 10:02:00 th Exp $


template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::open(const std::string &fname)
{
  output.open(fname.c_str());
  output.precision(14);
}

template<class Type, class StreamType>
SimpleDiagnostic<Type,StreamType>::~SimpleDiagnostic()
{
  output.close();
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


