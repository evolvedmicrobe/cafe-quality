using System;

namespace Bio
{
	public struct CigarElement
	{
		public int Length;
		public char Operation;

		public CigarElement (char operation, int length)
		{
			this.Length = length;
			this.Operation = operation;
		}
	}
}

