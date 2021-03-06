build_egg: 
	python setup.py bdist_egg

reinstall:
	-pip uninstall -y pygenesig 
	python setup.py sdist
	pip install --user dist/pygenesig-0.1.0.tar.gz

doc:
	cd docs/ && make html 

upload-doc: doc
	cd docs/_build/html && git add . && git add -u . && git commit -m "Update documentation. " && git push

