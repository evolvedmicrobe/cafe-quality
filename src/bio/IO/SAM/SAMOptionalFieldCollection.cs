using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using Bio.Util;
using System.IO;
using System.Globalization;

namespace Bio.IO.SAM
{
    /// <summary>
    /// This is a class that represents the optional field information
    /// 
    /// Only works on little endian architectures!!
    /// 
    /// </summary>
    public class SAMOptionalFieldCollection
    {
        /// <summary>
        /// A dictionary that gives the type of data and the location in the array where it starts.
        /// </summary>
        private Dictionary<string, KeyValuePair<char, int>> tagToDataTypeAndLocation;
        private byte[] data;
        
        /// <summary>
        /// Stores the metadata (read groups, etc). found at the end of BAM file.  Works 
        /// by directly copying the array and indexing all the read group information
        /// </summary>
        /// <param name="alignmentBlock">The original data from the parsed BAM</param>
        /// <param name="startIndex">where in the original data the tags start, this will be copied locally.</param>
        /// <param name="tagCache">A cache of strings to use for the tags, useful to avoid recreating a string from the two chars each time.</param>
        public SAMOptionalFieldCollection(byte[] alignmentBlock, int startIndex, Dictionary<ushort, string> tagCache = null)
        {
            tagToDataTypeAndLocation = new Dictionary<string, KeyValuePair<char, int>>(8);
            int sizeOfData = alignmentBlock.Length - startIndex;
            //now go through the array and for each tag, fill the data array with its value and the put the type in the dictionary
            //if (alignmentBlock.Length > startIndex + 4 && alignmentBlock[startIndex] != 0x0 && alignmentBlock[startIndex + 1] != 0x0)
            if (sizeOfData > 4 && alignmentBlock[startIndex] != 0x0 && alignmentBlock[startIndex + 1] != 0x0)
            {
                //copy all the tag data to a local array
                this.data = new byte[sizeOfData];
                Array.Copy(alignmentBlock, startIndex, data, 0, sizeOfData);
                //now index this array by the tag key and the position where the data starts        
                for (int origDataPosition = 0; origDataPosition < data.Length; )
                {
                    //read tag for block, getting cached string if possible
                    string tag;
                    if (tagCache != null)
                    {
                        ushort byte1 = data[origDataPosition];
                        ushort byte2 = data[origDataPosition + 1];
                        byte1 = (ushort)((byte1 << 8) + byte2);
                        bool cached = tagCache.TryGetValue(byte1, out tag);
                        if (!cached)
                        {
                            tag = System.Text.ASCIIEncoding.ASCII.GetString(data, origDataPosition, 2);
                            lock (tagCache) { tagCache[byte1] = tag; }
                        }
                    }
                    else
                    {
                        tag = System.Text.ASCIIEncoding.ASCII.GetString(data, origDataPosition, 2);
                    }
                    //now get the type of variable
                    origDataPosition += 2;
                    char vType = (char)data[origDataPosition++];
                    // SAM format supports [AifZH] for value type.
                    // In BAM, an integer may be stored as a signed 8-bit integer (c), unsigned 8-bit integer (C), signed short (s), unsigned
                    // short (S), signed 32-bit (i) or unsigned 32-bit integer (I), depending on the signed magnitude of the integer. However,
                    // in SAM, all types of integers are presented as type ʻiʼ. 
                    enterTagIntoIndex(vType, tag, ref origDataPosition);
                    //NOTE: Code previously here checked for valid value and threw an exception here, but this exception/validation is checked for in this method below, as while as when the value is set.
                }
            }
            else
            {
                //otherwise nothing
                this.data = new byte[0];
            }
        }

        /// <summary>
        /// Get the data for a tag
        /// </summary>
        /// <param name="key">The tag to get data for</param>
        /// <returns>The data or a </returns>
        public object this[string key]
        {
            get
            {
                if (key.Length != 2)
                {
                    throw new ArgumentException("All SAM Option fields must be of length 2, you asked for option field: " + key);
                }
                KeyValuePair<char, int> val;
                bool hasVal = this.tagToDataTypeAndLocation.TryGetValue(key, out val);
                if (!hasVal)
                    return null;
                else
                {
                    int statPos = val.Value;
                    return deSerializePosition(val.Key, ref statPos);
                }
            }
        }
        /// <summary>
        /// Convert this read only collection into a dictionary of option fields 
        /// </summary>
        /// <returns></returns>
        public Dictionary<string,SAMOptionalField> ConvertToDictionary()
        {              
                Dictionary<string,SAMOptionalField> toReturn=new Dictionary<string,SAMOptionalField>(tagToDataTypeAndLocation.Count);
                foreach(var kv in tagToDataTypeAndLocation)
                {
                    char type=kv.Value.Key;
                    object value=this[kv.Key].ToString();
                    string tag=kv.Key;
                    SAMOptionalField sof=new SAMOptionalField(tag,value,type);
                    toReturn[tag]=sof;
                }
                return toReturn;            
        }

        /// <summary>
        /// Get the data for a read tag from the underlying byte array.
        /// </summary>
        /// <param name="vType"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        private object deSerializePosition(char vType, ref int pos)
        {
            object obj;
            switch (vType)
            {
                case 'A':  //  Printable character
                    obj = (char)data[pos];
                    break;
                case 'c': //signed 8-bit integer                            
                    int value = (data[pos] & 0x7F);
                    if ((data[pos] & 0x80) == 0x80)
                    {
                        value= value + sbyte.MinValue;
                    }
                    obj = value;
                    break;
                case 'C'://uint8
                    obj = (int)data[pos];                    
                    break;
                case 's'://int16
                case 'S'://uint16
                    obj=BitConverter.ToInt16(data, pos);
                    break;
                case 'i'://int32
                case 'I'://uint32
                    obj = BitConverter.ToInt32(data, pos);
                    break;
                case 'f'://float
                    obj = BitConverter.ToSingle(data, pos);
                    break;
                case 'Z'://printable string 
                    int len = GetStringLength(data, pos);
                    obj = System.Text.ASCIIEncoding.ASCII.GetString(data, pos, len - 1);
                    break;
                case 'H'://byte array in hex format
                    len = GetStringLength(data, pos);
                    //obj = System.Text.ASCIIEncoding.ASCII.GetString(data, pos, len - 1);
                    obj = Helper.GetHexString(data, pos, len - 1);
                    break;
                case 'B'://integer or numeric array
                    char arrayType = (char)data[pos];
                    pos++;
                    int arrayLen = Helper.GetInt32(data, pos);
                    pos += 4;
                    StringBuilder strBuilder = new StringBuilder();
                    strBuilder.Append(arrayType);
                    for (int i = 0; i < arrayLen; i++)
                    {
                        strBuilder.Append(',');
                        string strValue = deSerializePosition(arrayType, ref pos).ToString();
                        strBuilder.Append(strValue);
                    }
                    obj = strBuilder.ToString();
                    break;
                default:
                    throw new FileFormatException(Properties.Resource.BAM_InvalidOptValType);
            }
            return obj;
        }    
   
        /// <summary>
        /// Index the tag elements in the data array and advance the index
        /// </summary>
        /// <param name="vType"></param>
        /// <param name="tag"></param>
        private void enterTagIntoIndex(char vType,string tag, ref int origDataPosition)
        {
            //some scenarios here, either we have to upcast to for bytes and copy
            //have to transfer four bytes, or have 
            int len;
                switch (vType)
                {
                    case 'A':  //  Printable character
                    case 'c': //signed 8-bit integer                    
                    case 'C'://uint8
                        tagToDataTypeAndLocation[tag] = new KeyValuePair<char, int>(vType, origDataPosition);    
                        origDataPosition++;
                        break;
                    case 's'://int16
                    case 'S'://uint16
                        tagToDataTypeAndLocation[tag] = new KeyValuePair<char, int>(vType, origDataPosition);    
                        origDataPosition+=2;
                        break;
                   case 'i'://int32
                   case 'I'://uint32
                   case 'f'://float
                        tagToDataTypeAndLocation[tag] = new KeyValuePair<char, int>(vType, origDataPosition);    
                        origDataPosition+=4;
                        break;
                   case 'Z'://printable string 
                   case 'H'://byte array in hex format
                        tagToDataTypeAndLocation[tag] = new KeyValuePair<char, int>(vType, origDataPosition);
                        len = GetStringLength(data, origDataPosition);
                        origDataPosition += len;
                        break;
                    case 'B'://integer or numeric array
                        tagToDataTypeAndLocation[tag] = new KeyValuePair<char, int>('B', origDataPosition);
                        char arrayType = (char)data[origDataPosition];
                        origDataPosition++;
                        int arrayLen = Helper.GetInt32(data, origDataPosition);
                        origDataPosition += 4;
                        origDataPosition += arrayLen;
                        break;
                    default:
                        throw new FileFormatException(Properties.Resource.BAM_InvalidOptValType);                
            }
        }

        /// <summary>
        /// Gets the length of the string in byte array.
        /// </summary>
        /// <param name="array">Byte array which contains string.</param>
        /// <param name="startIndex">Start index of array from which string is stored.</param>
        private int GetStringLength(byte[] array, int startIndex)
        {
            int i = startIndex;
            while (i < array.Length && array[i] != '\x0')
            {
                i++;
            }
            return i + 1 - startIndex;
        }

        /// <summary>
        /// Get the TAGs to append to the end of a BAM stream.
        /// </summary>
        /// <returns></returns>
        //public IReadOnlyCollection<byte> GetOriginalByteArray()
        //{
        //    return Array.AsReadOnly<byte>(data); 
        //}
    }
}
