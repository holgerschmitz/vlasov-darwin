#if defined INT
#define TYPE int
#define DF_TYPE DFNT_INT32
#elif defined FLOAT
#define TYPE float
#define DF_TYPE DFNT_FLOAT32
#elif defined DOUBLE
#define TYPE double
#define DF_TYPE DFNT_FLOAT64
#else
#error type unknown
#endif

template<int RANK>
HDFistream& HDFistream::operator>> (Matrix<TYPE, RANK>& m)
{
  if (sd_id < 0 || sds_cnt >= num_datasets)
    return *this;

  int32 rank, num_type, attr;
  NumArray<int32> dim_sizes(RANK);
  NumArray<int32> start(RANK); start.clear();
  char name[MAX_NC_NAME];
  
  int32 sds_id=SDselect(sd_id, sds_cnt++);

  SDgetinfo(sds_id, name, &rank, dim_sizes, &num_type, &attr);

  if (rank != RANK) {
    cerr << "HDFistream::operator>>: wrong dimension" << endl;
    return *this;
  }

  Index l(RANK), h(RANK);

  SDreadattr(sds_id,0,l);
  for (int d = 0; d < RANK; d++) {
    h[d] = l[d] + dim_sizes[RANK-d-1] - 1;
  }
  
  m.resize(l, h);

  if (num_type == DF_TYPE) {
    SDreaddata(sds_id, start, NULL, dim_sizes, m.data());
  } else if (num_type == DFNT_INT32) {
    int size = m.getSize();
    int *tmp = new int[size];
    SDreaddata(sds_id, start, NULL, dim_sizes, tmp);
    for (int i = 0; i < m.getSize(); i++) {
      m[i] = TYPE(tmp[i]);
    }
    delete[] tmp;
  } else if (num_type == DFNT_FLOAT32) {
    int size = m.getSize();
    float *tmp = new float[size];
    SDreaddata(sds_id, start, NULL, dim_sizes, tmp);
    for (int i = 0; i < m.getSize(); i++) {
      m[i] = TYPE(tmp[i]);
    }
    delete[] tmp;
  } else if (num_type == DFNT_FLOAT64) {
    int size = m.getSize();
    double *tmp = new double[size];
    SDreaddata(sds_id, start, NULL, dim_sizes, tmp);
    for (int i = 0; i < m.getSize(); i++) {
      m[i] = TYPE(tmp[i]);
    }
    delete[] tmp;
  } else {
    cerr << "HDFistream::operator>>: convert failed" << endl;
  }
  SDendaccess(sds_id);

  return *this;
}

template<int RANK>
HDFostream& HDFostream::operator<< (const Matrix<TYPE, RANK>& m)
{
  int32 sds_id;
  NumArray<int32> dim_sizes(RANK);
  for (int d = 0; d < RANK; d++)
    dim_sizes[d] = m.getDims(RANK-d-1);

  NumArray<int32> start(RANK); start.clear();

  if (flags & convert_to_int) {
    int size = m.getSize();
    int *tmp = new int[size];
    for (int i = 0; i < size; i++) {
      tmp[i] = int(m[i]);
    }
    sds_id = SDcreate(sd_id, "Matrix", DFNT_INT32, RANK, dim_sizes);
    SDsetattr(sds_id,"low",DFNT_INT32, RANK, (void *)m.getLow());
    SDsetdatastrs(sds_id, label, 0, 0, 0);
    SDwritedata(sds_id, start, 0, dim_sizes, tmp);
  } else if (flags & convert_to_float) {
    int size = m.getSize();
    float *tmp = new float[size];
    for (int i = 0; i < size; i++) {
      tmp[i] = float(m[i]);
    }
    sds_id = SDcreate(sd_id, "Matrix", DFNT_FLOAT32, RANK, dim_sizes);
    SDsetattr(sds_id,"low",DFNT_INT32, RANK, (void *)m.getLow());
    SDsetdatastrs(sds_id, label, 0, 0, 0);
    SDwritedata(sds_id, start, 0, dim_sizes, tmp);
    delete[] tmp;
  } else if (flags & convert_to_double) {
    int size = m.getSize();
    double *tmp = new double[size];
    for (int i = 0; i < size; i++) {
      tmp[i] = double(m[i]);
    }
    sds_id = SDcreate(sd_id, "Matrix", DFNT_FLOAT64, RANK, dim_sizes);
    SDsetattr(sds_id,"low",DFNT_INT32, RANK, (void *)m.getLow());
    SDsetdatastrs(sds_id, label, 0, 0, 0);
    SDwritedata(sds_id, start, 0, dim_sizes, tmp);
    delete[] tmp;
  } else {
    sds_id = SDcreate(sd_id, "Matrix", DF_TYPE, RANK, dim_sizes);
    SDsetattr(sds_id,"low",DFNT_INT32, RANK, (void *)m.getLow());
    SDsetdatastrs(sds_id, label, 0, 0, 0);
    SDwritedata(sds_id, start, 0, dim_sizes, (void *)m.data());
  }
  SDendaccess(sds_id);
  return *this;
}

#undef TYPE
#undef DF_TYPE
