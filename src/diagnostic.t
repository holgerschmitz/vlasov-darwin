// -*- C++ -*-
// $Id: $


template<class Type, class StreamType>
void SimpleDiagnostic<Type,StreamType>::open(const std::string &fname)
{
  output.open(fname.c_str());
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
