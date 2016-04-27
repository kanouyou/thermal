#ifndef IFdmException_HH
#define IFdmException_HH

#include <string>
#include <exception>

namespace FDM
{ class IFdmException; }

class FDM::IFdmException : public std::exception
{
  public:
    IFdmException(const std::string& aMessage);

    IFdmException(const std::string& aMessage, const char* aFile, const char* aFunction);

    virtual ~IFdmException() throw();

    const char* GetFileName() const;

    const char* GetFunction() const;

    virtual const char* what() const throw();    

  private:
    std::string fMessage;
    const char* fFileName;
    const char* fFunction;
};

#endif
