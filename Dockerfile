FROM python

# Fix timezone. Get user timezone from host with: $ timedatectl
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    llvm \
    make \
    tk-dev \
    wget \
    xz-utils \
    zlib1g-dev \
    seqtk \
    emboss \
    clustalw \
    ncbi-blast+ \
    git \
    python3-distutils \
    && rm -rf /var/lib/apt/lists/*
# Install python dependencies
RUN pip install --upgrade pip
RUN pip install --no-cache-dir biopython setuptools
# Clone the repository and CD into it
RUN git clone https://github.com/aditya-88/ASAP.git
WORKDIR /ASAP
# Link ASAP.py to the system path
RUN ln -s /ASAP/ASAP.py /usr/local/bin/ASAP
# Set permissions
RUN chmod +x /ASAP/ASAP.py
# Set the entrypoint
ENTRYPOINT ["ASAP"]

