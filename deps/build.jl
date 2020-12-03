using PyCall

if Sys.iswindows()
    @warn("""
                    FEniCS is not known to run on Windows!
                    Please try installing from Windows Subsystem for Linux.
          """)
end

try
    pyimport_conda("fenics", "fenics", "conda-forge")
    pyimport_conda("mshr", "mshr", "conda-forge")
catch ee
    typeof(ee) <: PyCall.PyError || rethrow(ee)
    @warn("""
                    Python Dependancies not installed
                    Please either:
                      - Rebuild PyCall using the path to FenicsPy using
                      - `ENV["PYTHON"]="... path of the python executable ..."; Pkg.build("PyCall"); Pkg.build("FenicsPy")`
          """)
end

