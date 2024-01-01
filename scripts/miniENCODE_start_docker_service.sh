#!/bin/bash

# miniENCODE Start Services Script
# Script Purpose: Initialize services for miniODP, including Apache, Shiny server, and SequenceServer.

# Obtain the host IP address
ip=$(hostname -I | awk '{print $1}')
if [ -z "$ip" ]; then
    echo "Failed to obtain IP address."
    exit 1
fi

# Update the HostIP in the miniODP configuration file
sed -i "s/\"HostIP\": \".*\"/\"HostIP\": \"http:\/\/$ip\"/g" /mnt/miniENCODE/miniODP/miniodp/static/SpeciesConfig.json

# Start Apache service
apache2ctl -D FOREGROUND > /dev/null 2>&1 &
sleep 1
if ! pgrep -x "apache2" > /dev/null;
then
    echo "Apache failed to start."
    exit 1
else
    echo "Apache started."
fi

# Start Shiny server
/usr/bin/shiny-server > /dev/null 2>&1 &
sleep 1
if ! pgrep -x "shiny-server" > /dev/null;
then
    echo "Shiny server failed to start."
    exit 1
else
    echo "Shiny server started."
fi

# Start SequenceServer
/usr/local/bin/sequenceserver -n 3 -p 4040 -d /db/blastdb > /dev/null 2>&1 &
sleep 1
if ! pgrep -x "sequenceserver" > /dev/null;
then
    echo "SequenceServer failed to start."
    exit 1
else
    echo "SequenceServer started."
fi

# Output the access information
echo -e "\nVisit the miniODP here: http://$ip/miniodp\n"
