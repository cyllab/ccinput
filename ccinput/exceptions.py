
# Base exception class for the project to make all exceptions easy to catch selectively
class CCInputException(Exception):
    pass

class InvalidParameter(CCInputException):
    pass

class InvalidXYZ(InvalidParameter):
    pass

class ImpossibleCalculation(CCInputException):
    pass

class MissingParameter(CCInputException):
    pass

class InternalError(CCInputException):
    pass

class UnimplementedError(CCInputException):
    pass
