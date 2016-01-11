/* 
 * File:   IndPos.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 26, 2015, 4:36 PM
 */

#ifndef INDPOS_H
#define	INDPOS_H

#include "Function.h"
#include "IndBox.h"


/**
 * \class IndPos
 * \brief %Indicator of the positive orthant
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on July 26, 2015, 4:36 PM
 * 
 * This class implements the indicator function of the positive orthant, that is
 * the cone
 *
 * \f[
 *    K_+ = \{x\in\mathbb{R}^n, x\geq 0\}
 * \f]
 *
 * The dimension \c n must be given to the constructor.
 * 
 * \ingroup Functions
 */
class IndPos : public IndBox {
public:
    virtual ~IndPos();
private:

};

#endif	/* INDPOS_H */

