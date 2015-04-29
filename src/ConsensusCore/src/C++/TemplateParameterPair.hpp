//
//  TemplateParameterPair.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/3/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__TemplateParameterPair__
#define __ConsensusCoreXCode__TemplateParameterPair__

#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <vector>

#include "TransitionParameters.hpp"
#include "ContextParameters.hpp"
#include "Mutation.hpp"
#include "Sequence.hpp"
#include "cassert"

namespace ConsensusCore {
    
    class WrappedTemplateParameterPair;
    
    class TemplateParameterPair {
        
    private:
        
        /* When a mutation takes place, it changes at most 2 other di-nucleotide contexts,
           the prior context and this one.  Rather than mutate the entire array holding the 
           template and transition parameters for a mutation we are just testing, we are going
           to just put the new contexts here, and change the function that acquires these 
           parameters for each.         
         */
        int mutantPosition;
        int mutationOffset; // 0, +1 or -1 for substitution, deletion, insertion.
        char mutantBP[2];
        TransitionParameters mutantParameters[2];
        
        void _ApplyMutationInPlace(const Mutation& mut, int start, const ContextParameters& ctx_params);
        
    public:
        std::string tpl;
        std::vector<TransitionParameters> trans_probs;
        
        // Constructors
        TemplateParameterPair(const std::string& tpl_, const std::vector<TransitionParameters>& trans_probs_);
        
        TemplateParameterPair();
        
        TemplateParameterPair(const std::string& tpl_, const ContextParameters& ctx);
        

        // Copy constructor
        // I want to delete this as it is a very expensive operation that I never want to happen
        // Reinstated for SWIG
        TemplateParameterPair(const TemplateParameterPair& other) = default;
        
        // Copy assignment - also deleting
        TemplateParameterPair& operator=(const TemplateParameterPair& rhs) = default;
        
        
        // Move assignment operator
        TemplateParameterPair& operator=(TemplateParameterPair&& rhs) = default;
        
        // Destructor
        ~TemplateParameterPair() = default;
        
        
        //Move Constructor
        TemplateParameterPair(TemplateParameterPair&& src) = default;
        
                TemplateParameterPair GetReverseComplement(const ContextParameters&ctx);
        
        /* Get a subsection of this template parameter pair */
        WrappedTemplateParameterPair GetSubSection(int start, int end);
        
        std::pair<char, TransitionParameters> GetTemplatePosition(int index) const {
            assert(index < tpl.length() && index >= 0);
            return std::pair<char, TransitionParameters>( tpl[index], trans_probs[index]);
        };
        
        // Has a virtual mutation been applied to the underlying template?
        bool VirtualMutationActive() const;
        
        /* Apply a "virtual" mutation to a template.  After applying, the bases and 
           and transition parameters for a position in the mutated template can be obtained with
           a call to GetVirtuallyMutatedTemplatePosition, but the underlying data for the original template is not
           copy or altered and can still be obtained by calling GetTemplatePosition
         */
        void ApplyVirtualMutation(const Mutation& mut, const ContextParameters& ctx_params);
        
        /*  Function gets the parameters for a position following a virtual mutation of a template.  A template can have a mutation
            "virtually" applied to it, and if so this will return the parameters from this virtually mutated template.
            (It is called virtual because rather than copying the underlying array and mutating it, we 
            simply adjust the array indexing to return the mutated values or original values.)
         */
        std::pair<char, TransitionParameters> GetVirtuallyMutatedTemplatePosition(int index) const {
            assert(index < tpl.length() - mutationOffset && index >= 0);
            // First we need to handle the case where the position overlaps with what has been mutated.
            if ( index < (mutantPosition - 1) ) {
                return GetTemplatePosition(index);
            }
            else if(index > mutantPosition) {
                return GetTemplatePosition(index + mutationOffset);
            }
            else {
                if (index == mutantPosition) {
                    auto toR = std::pair<char, TransitionParameters> (mutantBP[1], mutantParameters[1]);
                    return toR;
                }
                else {
                    return std::pair<char, TransitionParameters> (mutantBP[0], mutantParameters[0]);
                }
                
            }
            
        };
        
        /* 
            Get the length of a subsection of the template on
            the new "virtual" template, after accounting for any offsets
            introduced by the mutation.
         */
        int VirtualLength(int start, int length) const {
            int end = start + length; //exclusive
            if (mutantPosition >= start && mutantPosition < end) {
                return length - mutationOffset;
            }
            else {
                return length;
            }
        }
        
        
        void ClearVirtualMutation();
        
        void ApplyRealMutation(const Mutation& mut, const ContextParameters& ctx_params);
        void ApplyRealMutations(const std::vector<Mutation>& muts, const ContextParameters& ctx_params);

        
    };
    
    /*  Represents a class that wraps around a given template/parameter pair.
     Previously, each read had it's own template, which led to a nice design layout
     but a lot of duplication of data.  Now, there will only be one template for
     the forward and reverse positions, and this will be wrapped by these
     classes which will give access to the underlying data by converting the
     indexes directly.
     
     Wrapping and doing this index conversion is nasty, achtung.  It is however more performant */
    class WrappedTemplateParameterPair {
    
    private:
        TemplateParameterPair* base;
        int start;
        int length;
        
    public:
        
        int Length() const {
            return length;
        }
        
        int Start() const {
            return start;
        }
        
        
        int VirtualLength() const {
            return base->VirtualLength(start, length);
        }
        
        // Has a virtual mutation been applied to the underlying template?
        bool VirtualMutationActive() const;
        
        // Solely for the benefit of SWIG
        WrappedTemplateParameterPair() = default;
        
        WrappedTemplateParameterPair(const std::string& tpl_,
                                     const std::vector<TransitionParameters>& trans_probs_);
        
        WrappedTemplateParameterPair(TemplateParameterPair* base, int start, int length);
        
        // Copy constructor
        WrappedTemplateParameterPair(const WrappedTemplateParameterPair& other) = default;
        
        WrappedTemplateParameterPair(const std::string& tpl_, const ContextParameters& ctx);
        
        // Move assignment operator
        WrappedTemplateParameterPair& operator=(WrappedTemplateParameterPair&& rhs) = default;
        
        // Destructor
        ~WrappedTemplateParameterPair() = default;
        
        // Copy assignment
        WrappedTemplateParameterPair& operator=(const WrappedTemplateParameterPair& rhs) = default;
        
        //Move Constructor
        WrappedTemplateParameterPair(WrappedTemplateParameterPair&& src) = default;
        
        std::pair<char, TransitionParameters> GetTemplatePosition(int index) const {
            int newIndex = index + start;
            return base->GetTemplatePosition(newIndex);
        }
        
        std::pair<char, TransitionParameters> GetVirtuallyMutatedTemplatePosition(int index) const {
            int newIndex = index + start;
            return base->GetVirtuallyMutatedTemplatePosition(newIndex);
        }
        
    };

}

#endif /* defined(__ConsensusCoreXCode__TemplateParameterPair__) */
