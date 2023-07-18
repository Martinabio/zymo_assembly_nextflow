FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    git \
    pip

RUN pip3 install git+https://github.com/biowilko/Porechop.git

ENV PATH=/opt/bin:$PATH

CMD ["/bin/bash"]
