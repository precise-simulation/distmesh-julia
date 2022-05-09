rem Distmesh Benchmark Windows main test run script.

rem Copyright 2013-2018 Precise Simulation, Ltd.


SET OCTAVE="C:\Octave\Octave-4.0.3\bin\octave-cli.exe"
rem C:\Octave\Octave-4.2.1\octave.vbs
SET MATLAB="C:\Program Files\MATLAB\R2017a\bin\matlab.exe"
rem SET MATLAB="C:\Program Files\MATLAB\R2013b\bin\matlab.exe"
SET JULIA="%HOME%\AppData\Local\Julia-0.6.2\bin\julia.exe"


rm -f circ_*.txt
rm -f poly_*.txt
rm -f *.jpg


%OCTAVE% --no-gui distmesh_bench.m


%MATLAB% -nojvm -nosplash -r distmesh_bench


%JULIA% distmesh_bench.jl


%MATLAB% -nojvm -nosplash -r triangle_bench
%OCTAVE% --no-gui triangle_bench.m


%OCTAVE% --no-gui process_results.m
rem %MATLAB% -nosplash -r process_results


pause
