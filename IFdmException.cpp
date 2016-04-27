#include "IFdmException.h"

FDM::IFdmException::IFdmException(const std::string& aMessage)
    : fMessage(aMessage)
{}

FDM::IFdmException::IFdmException(const std::string& aMessage, const char* aFile,
                                  const char* aFunction)
    : fMessage(aMessage),
      fFileName(aFile),
      fFunction(aFunction)
{}

FDM::IFdmException::~IFdmException() throw()
{}

const char* FDM::IFdmException::what() const throw()
{ return fMessage.c_str(); }

const char* FDM::IFdmException::GetFileName() const
{ return fFileName; }

const char* FDM::IFdmException::GetFunction() const
{ return fFunction; }
