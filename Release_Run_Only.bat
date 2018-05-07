@echo off
cls

copy Output.txt Output_old.txt

"bin\Release\Simple_DE"
copy Output.txt Output_dispOnly.txt

FOR /L %%A IN (1,1,%1) DO (
	ECHO %%A
	"bin\Release\Simple_DE"
	copy Output.txt %2
)

echo(