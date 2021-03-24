# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TEST IMAGE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FROM blackfynn/base-processor-pandas-test:1-5800e81 as test

# for tests only
RUN cd /usr/bin && apk add --no-cache --update openjdk8 python2-dev gcc musl-dev libffi-dev jpeg-dev zlib-dev bash git openssh && \
    for i in $(find /usr/bin/../lib/jvm/default-jvm/bin -mindepth 1); do ln -s $i $(basename $i); done || true

RUN pip install --no-cache-dir javabridge==1.0.18 \
                               attrs==19.1.0 \
                               botocore==1.17.63 \
                               awscli==1.12.0 \
                               moto==1.3.7\
                               Pillow==5.3.0 \
                               pylint==1.9.3 \
                               pytest==4.0.1 \
                               pytest-cov==2.6.0 \
                               python-bioformats==1.5.2 \
                               boto3==1.7.84 \
                               pyaml==18.11.0

# code: processor
COPY base_image_microscopy_processor              /app/base_image_microscopy_processor
COPY brukertiff_processor/brukertiff_processor    /app/brukertiff_processor
COPY brukertiff_processor/run.py                  /app/run.py

# code: tests
COPY brukertiff_processor/tests                   /app/tests

WORKDIR /app

ENTRYPOINT [""]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRODUCTION IMAGE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FROM blackfynn/base-processor-pandas:1-5800e81 as prod

# for tests only
RUN cd /usr/bin && apk add --no-cache --update openjdk8 python2-dev gcc musl-dev libffi-dev openssl-dev jpeg-dev zlib-dev bash git openssh && \
    for i in $(find /usr/bin/../lib/jvm/default-jvm/bin -mindepth 1); do ln -s $i $(basename $i); done || true

RUN pip install --no-cache-dir javabridge==1.0.18 \
                               attrs==19.1.0 \
                               botocore==1.17.63 \
                               awscli==1.12.0 \
                               moto==1.3.7\
                               Pillow==5.3.0 \
                               pylint==1.9.3 \
                               pytest==4.0.1 \
                               pytest-cov==2.6.0 \
                               python-bioformats==1.5.2 \
                               boto3==1.7.84 \
                               pyaml==18.11.0

# code: processor
COPY base_image_microscopy_processor              /app/base_image_microscopy_processor
COPY brukertiff_processor/brukertiff_processor    /app/brukertiff_processor
COPY brukertiff_processor/run.py                  /app/run.py

# code: tests
COPY brukertiff_processor/tests                   /app/tests

WORKDIR /app

ENTRYPOINT [""]
