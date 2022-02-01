#ifndef SLIC3R_Platform_HPP
#define SLIC3R_Platform_HPP

<<<<<<< HEAD
=======
#include <string>

>>>>>>> master
namespace Slic3r {

enum class Platform
{
	Uninitialized,
	Unknown,
	Windows,
	OSX,
	Linux,
	BSDUnix,
};

enum class PlatformFlavor
{
	Uninitialized,
<<<<<<< HEAD
	Unknown,
	// For Windows and OSX, until we need to be more specific.
	Generic,
	// For Platform::Linux
	GenericLinux,
	LinuxOnChromium,
	// Microsoft's Windows on Linux (Linux kernel simulated on NTFS kernel)
	WSL,
	// Microsoft's Windows on Linux, version 2 (virtual machine)
	WSL2,
	// For Platform::BSDUnix
	OpenBSD,
=======
    Unknown,
    Generic,         // For Windows and OSX, until we need to be more specific.
    GenericLinux,    // For Platform::Linux
    LinuxOnChromium, // For Platform::Linux
    WSL,             // Microsoft's Windows on Linux (Linux kernel simulated on NTFS kernel)
    WSL2,            // Microsoft's Windows on Linux, version 2 (virtual machine)
    OpenBSD,         // For Platform::BSDUnix
    GenericOSX,      // For Platform::OSX
    OSXOnX86,        // For Apple's on Intel X86 CPU
    OSXOnArm,        // For Apple's on Arm CPU
>>>>>>> master
};

// To be called on program start-up.
void 			detect_platform();

Platform 		platform();
PlatformFlavor 	platform_flavor();

<<<<<<< HEAD
=======
std::string platform_to_string(Platform platform);
std::string platform_flavor_to_string(PlatformFlavor pf);

>>>>>>> master
} // namespace Slic3r

#endif // SLIC3R_Platform_HPP
