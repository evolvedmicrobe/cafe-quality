using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Bio.CrossPlatform
{
    public static class Environment    
	{
        public static bool RunningInMono
        {
            get { return Type.GetType("Mono.Runtime") != null; }
        }
		public enum Platform
		{
			Windows,
			Linux,
			Mac
		}
		public static Platform GetRunningPlatform()
		{
			switch (System.Environment.OSVersion.Platform)
			{
				case PlatformID.Unix:
				// Well, there are chances MacOSX is reported as Unix instead of MacOSX.
				// Instead of platform check, we'll do a feature checks (Mac specific root folders)
				if (Directory.Exists("/Applications")
				    & Directory.Exists("/System")
				    & Directory.Exists("/Users")
				    & Directory.Exists("/Volumes"))
					return Platform.Mac;
				else
					return Platform.Linux;
				case PlatformID.MacOSX:
					return Platform.Mac;
				default:
					return Platform.Windows;
			}
		}
		public static string GetFontDirectory()
		{
			string dir = "";
			Platform p = GetRunningPlatform ();
			if (p == Platform.Linux) {
				dir= GetLinuxFontDirectory ();
			} else if (p == Platform.Mac) {
				dir= "/Library/Fonts";
			}
			else if(p==Platform.Windows)
			{
				dir=System.Environment.GetEnvironmentVariable ("windir") + @"fonts";
			}
			if (!System.IO.Directory.Exists (dir)) {
				//throw new Exception ("Could not locate the font directory.  Tried to find fonts in " + dir + " assuming the platform was: " + p.ToString ());
			}
			return dir;

		}
        private static string GetLinuxFontDirectory()
        {
            string linuxDir="/usr/share/fonts/";
            if (RunningInMono && System.IO.Directory.Exists(linuxDir))
            {
                return linuxDir;
            }
            else
            {
                throw new Exception("Requested a linux font directory but the runtime is not mono.");
            }
        }

    }
}
