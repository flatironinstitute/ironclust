FROM matlabruntime/r2021b/release/update2/c0130000000000000

COPY ./irc2 /usr/bin/mlrtapp/

RUN apt-get update -y
RUN apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt-get install python3.8 -y
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 2
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

ENTRYPOINT ["/usr/bin/mlrtapp/irc2"]

