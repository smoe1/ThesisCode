#ifndef _EULERPARAMS_H_
#define _EULERPARAMS_H_

class IniDocument;
class EulerParams
{   
 public:
  double gamma;
  double x0;
  void init(IniDocument& ini_doc);
  void write_eulerhelp(const char* filename);
};
extern EulerParams eulerParams;

#endif
