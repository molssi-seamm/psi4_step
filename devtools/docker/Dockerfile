FROM molssi/mamba141

COPY ./psi4_step/data/seamm-psi4.yml /root/environment.yml

RUN mamba env update -f /root/environment.yml

WORKDIR /home
CMD ["psi4"]
