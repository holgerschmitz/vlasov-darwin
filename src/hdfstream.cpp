#include "hdfstream.h"

HDFstream::HDFstream()
{
  sd_id = -1; // no file opened;
  flags = 0;
  strcpy(label,"");
}

HDFstream::~HDFstream()
{
  if (sd_id >= 0)
    close();
}

void HDFstream::close()
{
  if (sd_id >= 0) {
    SDend(sd_id);
  }
  sd_id = -1;
}

int HDFstream::good() const
{
  if (sd_id < 0) {
    return 0;
  }
  return 1;
}

HDFstream::operator int () const
{
  return good();
}

int HDFstream::setf(int f)
{
  return flags |= f;
}

int HDFstream::unsetf(int f)
{
  return flags &= ~f;
}

void HDFstream::setLabel(const char* l)
{
  strncpy(label,l,256);
}

// ----------------------------------------------------------------------

HDFistream::HDFistream(const char* fname)
{
  open(fname);
}

int HDFistream::open(const char* fname)
{
  if (sd_id >= 0)
    close();
  if ((sd_id = SDstart(fname, DFACC_RDONLY)) < 0) {
    cerr << "HDFistream::open: " << fname << " failed" << endl;
    return 0;
  }
  sds_cnt = 0;
  int32 dummy;
  SDfileinfo(sd_id, &num_datasets, &dummy);
  return 1;
}

int HDFistream::select(int cnt)
{
  if (!HDFstream::good() || cnt < 0 || cnt >= num_datasets) {
    return 0;
  }
  sds_cnt = cnt;
  return 1;
}

int HDFistream::good() const
{
  if (!HDFstream::good() || sds_cnt < 0 || sds_cnt >= num_datasets) {
    return 0;
  }
  return 1;
}

int HDFistream::getNumDatasets() const
{
  return num_datasets;
}

void HDFistream::getDims(int32 &rank, int32 dims[]) const
{
  int32 sds_id = SDselect(sd_id, sds_cnt);
  int32 num_type, attr;
  char name[255];

  SDgetinfo(sds_id, name, &rank, dims, &num_type, &attr);
}

const char* HDFistream::getLabel() const
{
  int32 sds_id = SDselect(sd_id, sds_cnt);

  SDgetdatastrs(sds_id, (char*) label, 0, 0, 0, 255);
  return label;
}
  


// ----------------------------------------------------------------------

HDFostream::HDFostream(const char* fname)
{
  open(fname);
}

int HDFostream::open(const char* fname)
{
  if (sd_id >= 0)
    close();
  if ((sd_id = SDstart(fname, DFACC_CREATE)) < 0) {
    cerr << "HDFostream::open: " << fname << " failed" << endl;
    return 0;
  }

  return 1;
}

