include ../../config.h

pysander:
	@echo "[PYSANDER] setup.py install"
	$(VB)MSANDERHOME=$(BASEDIR) $(PYTHON) setup.py install $(PYTHON_INSTALL) \
           > setup.py.log

clean:
	/bin/rm -fr build/

skip:
	@echo ""
	@echo "Skipping installation of pysander."
	@echo ""

