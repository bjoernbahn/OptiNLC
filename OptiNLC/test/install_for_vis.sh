apt update
apt install -y python3-pip python3-venv
python3 -m venv /opt/venv
/opt/venv/bin/python -m pip install --upgrade pip
/opt/venv/bin/python -m pip install numpy matplotlib