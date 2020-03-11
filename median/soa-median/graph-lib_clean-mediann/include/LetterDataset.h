#ifndef __LETTER_DATASET_H__
#define __LETTER_DATASET_H__

#include "Dataset.h"
#include "LetterGraph.h"

class LetterDataset : public Dataset<CMUPoint, double, int>
{

protected:

  /**
   * Load a complete dataset, ie. a list of QAP Instances
   */
  void loadDS(const char* filename);

  /**
   * Load a QAP instance from a .dat file
   */
  void loadQAP(const char* filename);

public:

  LetterDataset(const char* filename);
  LetterDataset(): Dataset<CMUPoint, double, int>(){}

};


#endif //__LETTER_DATASET_H__
