# NUMERICAL METHODS FOR COMPUTER ENGINEERING
# @author: Miguel Blanco Godón

## Usage:
Input is provided through files. Each element of the matrix / vector must be separated by a blank character.
In case of error, an exception is raised and prompted in red colour.
This errors can appear because:
- Bad argument calling
- File not existing
- Bad type sent to the methods
- Some arithmetic exception
> Arithmetic exceptions, such as division by zero and math domain exceptions, may occur if the method data constraints are not met. For example, in Cholesky's factorization, 
> input matrix must be simmetric and definite positive).