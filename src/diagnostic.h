// -*- C++ -*-
// $Id$

#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include <iostream>
#include <list>
#include <string>
#include "stlpwrapper.h"
#include "task.h"

class DiagnosticInterface : public Rebuildable {
  private:
      std::string fname;
      std::string append;
      int interval;
      int t;
  public:
      DiagnosticInterface();
      void execute();
  protected:
      virtual void open(std::string &)=0;
      virtual void write()=0;
      virtual void close()=0;
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
  private:
      bool appending();
};

typedef PtrWrapper<DiagnosticInterface> pDiagnosticInterface;
typedef std::list<pDiagnosticInterface> DiagList;

class DiagnosticManager {
  private:
      static DiagnosticManager *theManager;
      DiagList diags;
  public:
      static DiagnosticManager& instance();
      void addDiagnostic(DiagnosticInterface*);
      void execute();
  private:
      DiagnosticManager();
      DiagnosticManager(DiagnosticManager&);
};

template<class Type, class StreamType>
class SimpleDiagnostic : public DiagnosticInterface {
  private:
      Type *field;
      StreamType output;
  public:
      void setField(Type*);
  protected:
      void open(std::string &);
      void write();
      void close();
};

#include "diagnostic.t"

#endif
