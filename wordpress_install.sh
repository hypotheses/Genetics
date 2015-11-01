#!/bin/sh
# 
# Wordpress Setup Script
# https://github.com/digitalocean/do_user_scripts/blob/master/Ubuntu-14.04/cms/wordpress.sh
# This script will install and configure Wordpress on
# an Ubuntu 14.04 droplet
export DEBIAN_FRONTEND=noninteractive;

# Generate root and wordpress mysql passwords
rootmysqlpass=`dd if=/dev/urandom bs=1 count=32 2>/dev/null | base64 -w 0 | rev | cut -b 2- | rev`
wpmysqlpass=`dd if=/dev/urandom bs=1 count=32 2>/dev/null | base64 -w 0 | rev | cut -b 2- | rev`
# Write passwords to file
echo "Root MySQL Password: $rootmysqlpass" > /root/passwords.txt;
echo "Wordpress MySQL Password: $wpmysqlpass" >> /root/passwords.txt;
# Update Ubuntu
apt-get update;
apt-get -y upgrade;
# Install Apache/MySQL
apt-get -y install apache2 php5 php5-mysql mysql-server mysql-client unzip;
# [BS] Install php5 -- may be the cause of bronken installation?
# [BS] apt-get -y install libapache2-mod-php5 php5-mcrypt
# Download and uncompress Wordpress
wget https://wordpress.org/latest.zip -O /tmp/wordpress.zip;
cd /tmp/;
unzip /tmp/wordpress.zip;
# Set up database user
/usr/bin/mysqladmin -u root -h localhost create wordpress;
/usr/bin/mysqladmin -u root -h localhost password $rootmysqlpass;
/usr/bin/mysql -uroot -p$rootmysqlpass -e "CREATE USER wordpress@localhost IDENTIFIED BY '"$wpmysqlpass"'";
/usr/bin/mysql -uroot -p$rootmysqlpass -e "GRANT ALL PRIVILEGES ON wordpress.* TO wordpress@localhost";

# Configure wordpress
cp /tmp/wordpress/wp-config-sample.php /tmp/wordpress/wp-config.php;
sed -i "s/'DB_NAME', 'database_name_here'/'DB_NAME', 'wordpress'/g" /tmp/wordpress/wp-config.php;
sed -i "s/'DB_USER', 'username_here'/'DB_USER', 'wordpress'/g" /tmp/wordpress/wp-config.php;
sed -i "s/'DB_PASSWORD', 'password_here'/'DB_PASSWORD', '$wpmysqlpass'/g" /tmp/wordpress/wp-config.php;

# original script install wordpress on /var/www/html
mkdir -p /var/www/html/wordpress
cp -Rf /tmp/wordpress/* /var/www/html/wordpress/.;
rm -f /var/www/html/index.html;
chown -Rf www-data:www-data /var/www/html; # essential to allow wordpress to download and install plug-in without ssh/ftp setup
a2enmod rewrite;
service apache2 restart;