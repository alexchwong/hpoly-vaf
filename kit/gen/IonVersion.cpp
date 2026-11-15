/* Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved */
#include "IonVersion.h"

using namespace std;

string IonVersion::GetMajor()     {return "5";}
string IonVersion::GetMinor()     {return "18";}
string IonVersion::GetRelease()   {return "6";}
string IonVersion::GetGitHash()   {return "5fe78dc2";}
string IonVersion::GetGitPath()   {return "local";}
string IonVersion::GetBuildNum()  {return "local";}
string IonVersion::GetBuildTag()  {return "local";}
string IonVersion::GetBuildHost() {return "local";}

string IonVersion::GetFullVersion(std::string module)
{
	string v = module;
	v += " (";
	v += "5";
	v += ".";
	v += "18";
	v += ".";
	v += "6";
	v += ")\n";
	v += "Build: ";
	v += "local";
	v += " (";
	v += "local";
	v += ")\nSource: ";
	v += "local";
	v += " (";
	v += "5fe78dc2";
	v += ")\n";
	return v;
}

string IonVersion::GetVersion()
{
	string v = "5";
	v += ".";
	v += "18";
   	return (v);
}



