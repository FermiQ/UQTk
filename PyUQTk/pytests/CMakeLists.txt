project (UQTk)

set(PYTHON_EXECUTABLE python)

SET( CMAKE_SWIG_OUTDIR "${PROJECT_BINARY_DIR}" )

# Add python tests and run without "make install"
configure_file( PyModTest.py "${CMAKE_SWIG_OUTDIR}/PyModTest.py" COPYONLY )
add_test( NAME PyModTest COMMAND ${PYTHON_EXECUTABLE} PyModTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyArrayTest.py "${CMAKE_SWIG_OUTDIR}/PyArrayTest.py" COPYONLY )
add_test( NAME PyArrayTest COMMAND ${PYTHON_EXECUTABLE} PyArrayTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyArrayTest2.py "${CMAKE_SWIG_OUTDIR}/PyArrayTest2.py" COPYONLY )
add_test( NAME PyArrayTest2 COMMAND ${PYTHON_EXECUTABLE} PyArrayTest2.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyQuadTest.py "${CMAKE_SWIG_OUTDIR}/PyQuadTest.py" COPYONLY )
add_test( NAME PyQuadTest COMMAND ${PYTHON_EXECUTABLE} PyQuadTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

#configure_file( PyBCSTest.py "${CMAKE_SWIG_OUTDIR}/PyBCSTest.py" COPYONLY )
#add_test( NAME PyBCSTest COMMAND ${PYTHON_EXECUTABLE} PyBCSTest.py WORKING_DIRECTORY $ #{CMAKE_SWIG_OUTDIR} )

configure_file( PyBCSTest1D.py "${CMAKE_SWIG_OUTDIR}/PyBCSTest1D.py" COPYONLY )
add_test( NAME PyBCSTest1D COMMAND ${PYTHON_EXECUTABLE} PyBCSTest1D.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyBCSTest2D.py "${CMAKE_SWIG_OUTDIR}/PyBCSTest2D.py" COPYONLY )
add_test( NAME PyBCSTest2D COMMAND ${PYTHON_EXECUTABLE} PyBCSTest2D.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyBADPTest.py "${CMAKE_SWIG_OUTDIR}/PyBADPTest.py" COPYONLY )
add_test( NAME PyBADPTest COMMAND ${PYTHON_EXECUTABLE} PyBADPTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyRegressionTest.py "${CMAKE_SWIG_OUTDIR}/PyRegressionTest.py" COPYONLY )
add_test( NAME PyRegressionTest COMMAND ${PYTHON_EXECUTABLE} PyRegressionTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )

configure_file( PyGalerkinTest.py "${CMAKE_SWIG_OUTDIR}/PyGalerkinTest.py" COPYONLY )
add_test( NAME PyGalerkinTest COMMAND ${PYTHON_EXECUTABLE} PyGalerkinTest.py WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR} )
