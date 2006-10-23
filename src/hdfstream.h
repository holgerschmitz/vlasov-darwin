#ifndef HDFSTREAM_H
#define HDFSTREAM_H
//-----------------------------------------------------------------------------
#include <hdf/mfhdf.h>
#include "matrix.h"
//-----------------------------------------------------------------------------
  /** @file hdfstream.h
	* @brief IO class for HDF
	*
	* IO classes for handling HDF files and streams
	*/
//-----------------------------------------------------------------------------
//HDFstream
/** @brief IO class for handling HDF files
  *  
  * This is the abstract base class for HDF-IO- classes.
  * Implements the basic operations on HDFstreams as virtual methods.
  */
class HDFstream {
  protected:
  	/// ID
      int32 sd_id;
      ///Flags
      int flags;
      ///label
      char label[256];

  public:
	/// flags
      enum {
		convert_to_int    = 1, ///< convert to to integer
		convert_to_float  = 2, ///< convert to flaoting point
		convert_to_double = 4  ///< convert to doubel precision
      };
	/// constructor 
      HDFstream();
	///copy constructur 
      HDFstream(const HDFstream&);
	/// destructor
      virtual ~HDFstream();

	/// open  file 
      virtual int open(const char*)=0;
	/// close file 
      virtual void close();
	/// return true=1 if data are still available 
      virtual int good() const;
      /// return true=1 if data are still available 
      operator int () const;
	/// assign 
      HDFstream& operator = (const HDFstream&);
	/// set flags 
      int setf(int);
	/// unset flags 
      int unsetf(int);
      /// add label 
      void setLabel(const char *label);
};
//HDFstream
//---------------------------------------------------------------
//HDFistream
/** @brief Input stream for HDF files */
class HDFistream : public HDFstream {
	int32 sds_cnt,       ///< counter for the send sets
	      num_datasets;  ///< total number of datasets

  public:
	/// constructor 
      HDFistream() : HDFstream() {}
	/// copy constructor */
      HDFistream(const HDFstream&);
      /// constructor, opens HDF file "fname", selects first dataset 
      HDFistream(const char* fname);

  /// opens HDF file "fname", selects first dataset 
  virtual int open(const char*);
  /// return true=1 if data are still available 
  virtual int good() const;
  /// selects HDF dataset with number "cnt" 
  int select(int cnt);
  /// return number of datasets 
  int getNumDatasets() const;
  /// return rank and dims of currently selected field 
  void getDims(int32 &rank, int32 dims[]) const;
  /// get label of currently selected field 
  const char* getLabel() const;

  /// stream input operator for a matrix of int's
  template<int RANK>
  HDFistream& operator>> (Matrix<int, RANK>&);
  /// stream input operator for a matrix of float's
  template<int RANK>
  HDFistream& operator>> (Matrix<float, RANK>&);
  /// stream input operator for a matrix of doubles's
  template<int RANK>
  HDFistream& operator>> (Matrix<double, RANK>&);
};
//HDFistream
//-----------------------------------------------------------------------------
//HDFostream
/** @brief output stream for HDF files */
class HDFostream : public HDFstream {
  public:
	/// constructor 
      HDFostream() : HDFstream() {}
      /// copy constructor 
      HDFostream(const HDFstream&);
      /// constructor, opens HDF file "fname" 
      HDFostream(const char* fname);

	/// open file 
      virtual int open(const char*);
	 /// set floating point/double precision, does nothing (compablility?)
      void precision(int) {}
      
	/// stream output operator for a matrix of int's
      template<int RANK>
      HDFostream& operator<< (const Matrix<int, RANK>&);
	/// stream output operator for a matrix of float's
      template<int RANK>
      HDFostream& operator<< (const Matrix<float, RANK>&);
	/// stream output operator for a matrix of double's
      template<int RANK>
      HDFostream& operator<< (const Matrix<double, RANK>&);
};
//HDFostream
//-----------------------------------------------------------------------------
///include the implementation for int type, outdated hack, not unlike template instancation
#define INT
#include "hdfstream.t"
#undef INT
///include the implementation for float type, outdated hack, not unlike template instancation
#define FLOAT
#include "hdfstream.t"
#undef FLOAT
///include the implementation for double type, outdated hack, not unlike template instancation
#define DOUBLE
#include "hdfstream.t"
#undef DOUBLE
//-----------------------------------------------------------------------------
#endif // HDFSTREAM_H
