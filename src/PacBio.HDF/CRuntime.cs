#region Copyright (c) 2010, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// THIS SOFTWARE CONSTITUTES AND EMBODIES PACIFIC BIOSCIENCES CONFIDENTIAL
// AND PROPRIETARY INFORMATION.
//
// Disclosure, redistribution and use of this software is subject to the
// terms and conditions of the applicable written agreement(s) between you
// and Pacific Biosciences, where you refers to you or your company or
// organization, as applicable.  Any other disclosure, redistribution or
// use is prohibited.
//
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#endregion

using System;
using System.Runtime.InteropServices;
using System.IO;

namespace PacBio.HDF
{
    /// <summary>
    /// DLL imports for low level c runtime calls
    /// </summary>
    /// When using C libraries sometimes we need to call things like 'free' to keep from leaking memory in the C heap
    public class CRuntime
    {
        // Platform check code from
        // http://stackoverflow.com/questions/10138040/how-to-detect-properly-windows-linux-mac-operating-systems

        public enum Platform
        {
            Windows,
            Linux,
            Mac
        }

        private static Platform GetRunningPlatform()
        {
            switch (Environment.OSVersion.Platform)
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

        public static readonly Platform RunningPlatform = GetRunningPlatform();


        //The actual DLL this comes from is msvcrt80 but it seems that the .net loader knows this 'magic' name
        
        private const string LinuxRuntimeName = "libc.so.6";

        /// <summary>
        /// C library free() routine
        /// </summary>
        /// <param name="addr"></param>
        /// <returns></returns>
        [DllImport(LinuxRuntimeName, EntryPoint = "free", CallingConvention = CallingConvention.Cdecl)]
        public static extern void free_Linux(IntPtr addr);


        private const string WindowsRuntimeName = "hdf5.dll";

        /// <summary>
        /// C library free() routine. Same heap as the one used by hdf5
        /// </summary>
        /// <param name="addr"></param>
        /// <returns></returns>
        [DllImport("msvcr100.dll", EntryPoint = "free", CallingConvention = CallingConvention.Cdecl)]
        public static extern void free_Windows(IntPtr addr);


        private const string MacRuntimeName = "libc.dylib";

        [DllImport(MacRuntimeName, EntryPoint = "free", CallingConvention = CallingConvention.Cdecl)]
        public static extern void free_Mac(IntPtr addr);


        public static void free(IntPtr addr)
        {
            if (RunningPlatform == Platform.Windows)
            {
                free_Windows(addr);
            }
            else if (RunningPlatform == Platform.Linux)
            {
                free_Linux(addr);
            }
            else
            {
                free_Mac(addr);
            }
        }

        /// <summary>
        /// dlopen in linux
        /// </summary>
        /// <param name="dllToLoad"></param>
        /// <returns></returns>
        public static IntPtr GetModuleHandle(string dllToLoad)
        {
            IntPtr hdl = IntPtr.Zero;
            if (RunningPlatform == Platform.Windows)
            {
                hdl = GetModuleHandleInternal_Win(dllToLoad);
            }
            else if (RunningPlatform == Platform.Linux)
            {
                // Console.WriteLine((new System.Diagnostics.StackTrace(true)).ToString());

                // Note -- flags==1 is RTLD_LAZY
                // flags 0x100 is RTLD_GLOBAL
                hdl = GetModuleHandleInternal_Linux(dllToLoad, 1 | 0x100);
                
                if (hdl == IntPtr.Zero)
                    throw new Exception("DLL load failure on " + dllToLoad + ": " + GetDLLLoadError_Linux());
            }
            else // Mac 
            {
                //                #define RTLD_LAZY       0x1
                // ..
                //                #define RTLD_GLOBAL     0x8
                hdl = GetModuleHandleInternal_Mac(dllToLoad, 1 | 0x8);

                if (hdl == IntPtr.Zero)
                    throw new Exception("DLL load failure on " + dllToLoad + ": " + GetDLLLoadError_Mac());

            }

            if (hdl == IntPtr.Zero)
                throw new Exception("Failed to GetModuleHandle for: " + dllToLoad);

            return hdl;
        }

        [DllImport("kernel32.dll", EntryPoint = "GetModuleHandle")]
        private static extern IntPtr GetModuleHandleInternal_Win(string dllToLoad);

        [DllImport("libdl.so.2", EntryPoint = "dlopen")]
        private static extern IntPtr GetModuleHandleInternal_Linux(string dllToLoad, int flags);

        [DllImport("libdl.so.2", EntryPoint = "dlerror")]
        private static extern string GetDLLLoadError_Linux();

        [DllImport("libSystem.B.dll", EntryPoint = "dlopen")]
        private static extern IntPtr GetModuleHandleInternal_Mac(string dllToLoad, int flags);

        [DllImport("libSystem.B.dll", EntryPoint = "dlerror")]
        private static extern string GetDLLLoadError_Mac();


		
        /// <summary>
        /// dlsym in linux
        /// </summary>
        /// <param name="hModule"></param>
        /// <param name="procedureName"></param>
        /// <returns></returns>
        public static IntPtr GetProcAddress(IntPtr hModule, string procedureName)
        {
            if (RunningPlatform == Platform.Windows)
                return GetProcAddressInternal_Win(hModule, procedureName);
            else if (RunningPlatform == Platform.Linux)
                return GetProcAddressInternal_Linux(hModule, procedureName);
            else
                return GetProcAddressInternal_Mac(hModule, procedureName);
        }

        [DllImport("kernel32.dll", EntryPoint = "GetProcAddress")]
        private static extern IntPtr GetProcAddressInternal_Win(IntPtr hModule, string procedureName);

        [DllImport("libdl.so.2", EntryPoint = "dlsym")]
        private static extern IntPtr GetProcAddressInternal_Linux(IntPtr hModule, string procedureName);

        [DllImport("libSystem.B.dylib", EntryPoint = "dlsym")]
        private static extern IntPtr GetProcAddressInternal_Mac(IntPtr hModule, string procedureName);


		
        /// <summary>
        /// dlclose in linux
        /// </summary>
        /// <param name="hModule"></param>
        /// <returns></returns>
        [DllImport("kernel32.dll")]
        public static extern bool FreeLibrary(IntPtr hModule);
    }
}
