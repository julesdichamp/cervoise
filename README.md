# CERVOISE

CERVOISE = Covid19 : Estimation of Regional Virus Organization Infered using Stochastic Equation

Notes can be found here: https://onedrive.live.com/edit.aspx?cid=cd3df43e61c65f01&page=view&resid=CD3DF43E61C65F01!514&parId=CD3DF43E61C65F01!513&app=Word

Contains:
- README.md
- Meetings.txt with past meetings notes
- config.ini file with parameters for initial configuration
- data directory
- cervoise.py (main)
- getFranceGeomDpts.py: compute data between and inside departments
- models.py: compute models 
- plotModels.py: plot models results

Install:
 - install pip (typically : *conda install pip*)
 - REQUIRES_PYTHON = '>=3.6.0'
 - from the cervoise folder : *python setup.py install*

TODO:
- Get data from cases/deaths 
- Get connections between regions (surface of contact) 
- Create graph (including those datas) 
- Chemical master equation + "diffusion" master equation
- Solver 
- Fit missing parameters 
- Validation procedure (either by training on Italy data) or by iteratively fit on n data and predict on n+1 
