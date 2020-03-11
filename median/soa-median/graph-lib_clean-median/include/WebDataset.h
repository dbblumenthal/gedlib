#ifndef __Web_DATASET_H__
#define __Web_DATASET_H__

#include "Dataset.h"
#include "WebGraph.h"

class WebDataset : public Dataset<WebNAtt, WebEAtt, int>
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

  WebDataset(const char* filename);
  WebDataset(): Dataset<WebNAtt, WebEAtt, int>(){}

};


#endif //__Web_DATASET_H__
