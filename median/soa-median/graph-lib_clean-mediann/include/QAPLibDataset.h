#ifndef __QAPLIB_DATASET_H__
#define __QAPLIB_DATASET_H__

#include "Dataset.h"
#include "QAPLibGraph.h"


class QAPLibDataset : public Dataset <int, int, int>
{

private:

  /**
   * @param filename  Dataset .ds
   */
  void loadDS(const char* filename);

  /**
   * @param filename  Instance of QAPLib
   */
  void loadQAP(const char* filename);

public:

  QAPLibDataset(const char* filname);

};

#endif //  __QAPLIB_DATASET_H__
