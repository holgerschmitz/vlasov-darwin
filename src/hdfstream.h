#ifndef HDFSTREAM_H
#define HDFSTREAM_H

#include <hdf/mfhdf.h>
#include "matrix.h"

/** IO for HDF files */
class HDFstream {
protected:
  int32 sd_id;
  int flags;
  char label[256];

public:
  /** flags*/
  enum {
    /** */
    convert_to_int    = 1,
    /** */
    convert_to_float  = 2,
    /** */
    convert_to_double = 4
  };

  /** constructor */
  HDFstream();

  /** copy constructur */
  HDFstream(const HDFstream&);

  /** destructor */
  virtual ~HDFstream();

  /** */
  virtual int open(const char*)=0;
  /** */
  virtual void close();
  /** return true=1 if data are still available */
  virtual int good() const;
  /** return true=1 if data are still available */
  operator int () const;
  /** assign */
  HDFstream& operator = (const HDFstream&);
  /** set flags */
  int setf(int);
  /** unset flags */
  int unsetf(int);
  /** add label */
  void setLabel(const char *label);
};

/** input stream for HDF files */
class HDFistream : public HDFstream {
  int32 sds_cnt, num_datasets;

public:
  /** constructor */
  HDFistream() : HDFstream() {}
  /** copy constructor */
  HDFistream(const HDFstream&);
  /** constructor, opens HDF file "fname", selects first dataset */
  HDFistream(const char* fname);

  /** opens HDF file "fname", selects first dataset */
  virtual int open(const char*);
  /** return true=1 if data are still available */
  virtual int good() const;
  /** selects HDF dataset with number "cnt" */
  int select(int cnt);
  /** return number of datasets */
  int getNumDatasets() const;
  /** return rank and dims of currently selected field */
  void getDims(int32 &rank, int32 dims[]) const;
  /** get label of currently selected field */
  const char* getLabel() const;

  /** */
  template<int RANK>
  HDFistream& operator>> (Matrix<int, RANK>&);
  /** */
  template<int RANK>
  HDFistream& operator>> (Matrix<float, RANK>&);
  /** */
  template<int RANK>
  HDFistream& operator>> (Matrix<double, RANK>&);
};

/** output stream for HDF files */
class HDFostream : public HDFstream {
public:
  /** constructor */
  HDFostream() : HDFstream() {}
  /** copy constructor */
  HDFostream(const HDFstream&);
  /** constructor, opens HDF file "fname" */
  HDFostream(const char* fname);

  /** */
  virtual int open(const char*);

  /** */
  template<int RANK>
  HDFostream& operator<< (const Matrix<int, RANK>&);
  /** */
  template<int RANK>
  HDFostream& operator<< (const Matrix<float, RANK>&);
  /** */
  template<int RANK>
  HDFostream& operator<< (const Matrix<double, RANK>&);
};

#define INT
#include "hdfstream.t"
#undef INT
#define FLOAT
#include "hdfstream.t"
#undef FLOAT
#define DOUBLE
#include "hdfstream.t"
#undef DOUBLE

#endif
