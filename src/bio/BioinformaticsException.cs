using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio
{
    public class BioinformaticsException : Exception
    {
        public BioinformaticsException(string message, Exception inner)
            : base(message, inner)
        { }
        public BioinformaticsException(string message) : base(message) { }

    }
}
