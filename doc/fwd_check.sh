#! /bin/sh
# $Id: fwd_check.sh,v 1.75 2011/07/20 16:04:48 lavr Exp $
# Author:   Denis Vakatov (vakatov@ncbi,nlm.nih.gov)
# Modified: Anton Lavrentiev (lavr@ncbi.nlm.nih.gov)
#
# Check for the status of FWDAEMON on the dispatching NCBI servers

delay_sec="$1"
delay_sec=${delay_sec:="10"}
netcat="`which netcat 2>/dev/null`"
temp="/tmp/`basename $0`.$$.tmp"
helper="./fwd_failure_helper.exe"
test -z "$netcat"  &&  netcat="`whereis netcat | sed 's/^[^:]*://;s/ //g'`"

cat <<EOF
http://www.ncbi.nlm.nih.gov/IEB/ToolBox/NETWORK/firewall.html

Checking connections to NCBI Firewall Daemons as of `date -u +'%b %d %Y %R GMT'`:
EOF

{
cat <<EOF
;130.14.25.13	5555	RETIRED
10.10.150.44	5555	INTERNAL
130.14.24.219   5555    INTERNAL
130.14.29.112	5860	RESERVED
130.14.29.112	5861	OK
130.14.29.112	5862	RESERVED
130.14.29.112	5863	RESERVED
130.14.29.112	5864	RESERVED
130.14.29.112	5865	RESERVED
130.14.29.112	5866	RESERVED
130.14.29.112	5867	OK
130.14.29.112	5868	OK
130.14.29.112	5869	RESERVED
130.14.29.112	5870	OK
130.14.29.112	4444	RETIRED
130.14.29.112	4445	OK
130.14.29.112	4446	RESERVED
130.14.29.112	4447	RESERVED
130.14.29.112	4448	RESERVED
130.14.29.112	4449	RESERVED
130.14.29.112	4450	RESERVED
130.14.29.112	4451	OK
130.14.29.112	4452	OK
130.14.29.112	4453	RESERVED
130.14.29.112	4454	OK
130.14.29.112	443	FB-OK
130.14.29.112	22	FB-OK
165.112.7.12	5860	RESERVED
165.112.7.12	5861	RESERVED
165.112.7.12	5862	RESERVED
165.112.7.12	5863	RESERVED
165.112.7.12	5864	OK
165.112.7.12	5865	OK
165.112.7.12	5866	OK
165.112.7.12	5867	RESERVED
165.112.7.12	5868	RESERVED
165.112.7.12	5869	RESERVED
165.112.7.12	5870	RESERVED
165.112.7.12	4444	RETIRED
165.112.7.12	4445	RESERVED
165.112.7.12	4446	RESERVED
165.112.7.12	4447	RESERVED
165.112.7.12	4448	OK
165.112.7.12	4449	OK
165.112.7.12	4450	OK
165.112.7.12	4451	RESERVED
165.112.7.12	4452	RESERVED
165.112.7.12	4453	RESERVED
165.112.7.12	4454	RESERVED
EOF
} |
while read x_host x_port x_status ; do
  test "`echo $x_host | grep -c '^[;]'`" != "0"  &&  continue
  if [ "$x_port" -lt "5860" -o "$x_port" -gt "5870" ]; then
    if [ "$x_port" -ne "22" -a "$x_port" -ne "443" ]; then
      test "$x_status" = "RESERVED"  &&  continue
      if [ _"$HTTP_CAF" = _"" -o _"$HTTP_CAF_EXTERNAL" != _"" ]; then
        test _"$HTTP_NCBI_RELAY" = _"" &&  continue
      fi
    fi
  fi
  if [ "$x_status" = "RETIRED" -o "$x_status" = "RESERVED" -o "$x_status" = "PENDING" ]; then
    printf '%s\t%-8s\n'		"${x_host}:${x_port}"	"$x_status"
    continue
  fi
  unset fb_port
  test "$x_status" = "FB-OK" && fb_port="TRUE"
  test "$x_status" = "READYING"  &&  unset x_status
  ( echo ; test -z "$netcat" && sleep ${delay_sec} ) | ${netcat:-telnet} $x_host $x_port >$temp 2>&1 &
  pid=$!
  trap 'rm -f $temp; kill $pid >/dev/null 2>&1' 1 2 15
  ( sleep `expr $delay_sec + 1`  &&  kill $pid ) >/dev/null 2>&1 &
  guard=$!
  wait $pid >/dev/null 2>&1
  kill $guard >/dev/null 2>&1
  test _"$HTTP_CAF_EXTERNAL" != _""  ||  \
    cp="`tail +3 $temp 2>/dev/null | grep -s '[[]\{0,1\}[0-9]\{1,3\}[.][0-9]\{1,3\}[.][0-9]\{1,3\}[.][0-9]\{1,3\}[]]\{0,1\}:[0-9]\{1,5\}'`"
  grep -qs 'NCBI Firewall Daemon:  Invalid ticket\.  *Connection closed\.' $temp >/dev/null 2>&1
  if   [ $? -eq 0 ]; then
    printf '%s\t%-8s%s\n'	"${x_host}:${x_port}"	"${x_status:-OKAY}"	"${cp:+	}${cp}"
  elif [ -z "$x_status" ]; then
    printf '%s\t%-8s\n'		"${x_host}:${x_port}"	"READYING"
  else
    x_reason="( telnet $x_host $x_port )"
    if [ -x "$helper" ]; then
      x_status="`$helper $x_host $x_port`"
      test $? -eq 0  &&  x_reason=''
    else
      x_status=''
    fi
    status=${x_status:-FAILED}
    test ! -z "$fb_port" && status="FB-$status"
    printf '%s\t%-8s%s\n'	"${x_host}:${x_port}"	"${status}"	"${x_reason:+	}${x_reason}"
  fi
  rm -f $temp
done 2>&1 | grep -v 'Terminated'
