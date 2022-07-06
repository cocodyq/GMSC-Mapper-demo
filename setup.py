from setuptools import setup, find_packages

setup(
    name="GMSC-mapper",  
    version="0.0.1_alpha", 
    description="Annotate smORFs based on the Global Microbial smORFs Catalog (GMSC)",  
    long_description=open("./README.md", "r").read(),  
    long_description_content_type="text/markdown",  
    url="https://github.com/cocodyq/GMSC-Mapper-demo",  
    author="Yiqian Duan",
    author_email="yqduan20@fudan.edu.cn",  
    license="MIT"
    classifiers=[  
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
    ],
    package_dir={"GMSC-Mapper-demo": "GMSC-mapper"}, 
    packages=[GMSC-mapper,GMSC-mapper.data],  
	include_package_data=True,
    python_requires=">=3.8",
    install_requires=open("./requirements.txt", "r").read().splitlines(),
    package_data={"GMSC-Mapper-demo": ["GMSC-mapper.data"]},
    zip_safe=False,
    entry_points={  
        "console_scripts": [
            "GMSC-mapper=GMSC-mapper.cli:main",
        ],
    },
)
