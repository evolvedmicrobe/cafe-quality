using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public class BaseAndQualityAndPosition
    {
        public int Position;
        
		public byte Base 
        { 
            get{return BaseWithQuality.Base;}
        }
        
		public ushort InsertionOffSet;

        public bool PositionIsInsertion
        {
            get 
            {
                return InsertionOffSet > 0;            
            }
        }         

        public BaseAndQuality BaseWithQuality;
        
		public BaseAndQualityAndPosition(int position, int insertSizeRelativeToPosition, BaseAndQuality BandQ)
        {
            this.BaseWithQuality = BandQ;
            this.Position = position;
            this.InsertionOffSet = (ushort) insertSizeRelativeToPosition;
        }
    }
}
