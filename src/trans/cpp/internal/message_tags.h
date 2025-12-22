// (C) Copyright 2025- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#ifndef MESSAGE_TAGS_H
#define MESSAGE_TAGS_H

// Tags for identifying MPI messages
namespace message_tags {
    const int mtagletr   = 18000;
    const int mtagletr   = 18000;
    const int mtagml     = 19000;
    const int mtaglg     = 20000;
    const int mtagpart   = 21000;
    const int mtagdistsp = 22000;
    const int mtaggl     = 23000;
    const int mtaglm     = 24000;
    const int mtagdistgp = 25000;
}

#endif // MESSAGE_TAGS_H
