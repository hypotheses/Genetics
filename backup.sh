#!/bin/bash -w

# Running with sudo 
# To automate the back up of a web server 
# 1) on remote backup storage server
#--1 create `backup' user
# the home of this user is /var/backups
#--2 create a folder backups 

# crontab file -- Uncomment line 30-34 
# Edit this file to introduce tasks to be run by cron.
#
# Each task to run has to be defined through a single line
# indicating with different fields when the task will be run
# and what command to run for the task
#
# To define the time you can provide concrete values for
# minute (m), hour (h), day of month (dom), month (mon),
# and day of week (dow) or use '*' in these fields (for 'any').#
# Notice that tasks will be started based on the cron's system
# daemon's notion of time and timezones.
#
# Output of the crontab jobs (including errors) is sent through
# email to the user the crontab file belongs to (unless redirected).
#
# For example, you can run a backup of all your user accounts
# at 5 a.m every week with:
# 0 5 * * 1 tar -zcf /var/backups/home.tgz /home/
#
# For more information see the manual pages of crontab(5) and cron(8)
#
# m h  dom mon dow   command
# ---------------------------------
MAILTO=bhoom.suk@mahidol.ac.th
30 3 * * * /bin/tar -czf /backup/www_initial_backup-`date '+%m%d%y'`.tar.gz /var/www; scp -i /home/bhoom/.ssh/id_rsa_backup /backup/www_initial_backup-`date '+%m%d%y'`.tar.gz backup@10.90.202.239:/var/backups/backups
0 4 * * * mysqldump --databases mysql redcap test --events -uroot -predcapdatastorage | gzip > /backup/mysql/sql-01--`date '+%m%d%y'`.sql.gz \
/usr/bin/scp -i /home/bhoom/.ssh/id_rsa_backup /backup/mysql/sql-01-`date +\%m\%d\%y`.sql.gz  backup@10.90.202.239:~/backups
5 4 * * * mysqldump --single-transaction information_schema -uroot -predcapdatastorage | gzip > /backup/mysql/sql-02--`date '+%m%d%y'`.sql.gz \
/usr/bin/scp -i /home/bhoom/.ssh/id_rsa_backup /backup/mysql/sql-02-`date +\%m\%d\%y`.sql.gz  backup@10.90.202.239:~/backups
10 4 * * * mysqldump --single-transaction performance_schema -uroot -predcapdatastorage | gzip > /backup/mysql/sql-03--`date '+%m%d%y'`.sql.gz \
/usr/bin/scp -i /home/bhoom/.ssh/id_rsa_backup /backup/mysql/sql-03-`date +\%m\%d\%y`.sql.gz  backup@10.90.202.239:~/backups

