// -*- C++ -*-
// $Id$

#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include <iostream>
#include <list>
#include <string>
#include "stlpwrapper.h"
#include "task.h"

/** Interface for diagnostic tasks.
 * 
 * This interface can be used to implement different types of diagnostics.
 * The DiagnosticInterface is closely related to the DiagnosticManager. When an
 * instance of the interface is created it will register itself with the 
 * DiagnosticManager. This then takes the responsibility of calling the execute
 * method of the DiagnosticInterface.
 */
class DiagnosticInterface : public Rebuildable {
  private:
      /// The file name into which to write
      std::string fname;
      /** A string specifying whether to append to the file or whether to
       * write into a new file at each turn.
       * Only the first letter is checked for equality to 'y'
       */
      std::string append;
      int interval;
      int t;
  public:
      DiagnosticInterface();
      void execute();
  protected:
      virtual void open(const std::string &)=0;
      virtual void write()=0;
      virtual void close()=0;
      virtual bool singleOut() { return false; }
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
  private:
      bool appending();
      std::string parsedFileName();
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
      SimpleDiagnostic() {}
      void setField(Type*);
  protected:
      void open(const std::string &);
      void write();
      void close();
};


#include "diagnostic.t"

#endif
